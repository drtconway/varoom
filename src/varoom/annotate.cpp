#include "varoom/command.hpp"

#include <chrono>
#include <map>
#include <iostream>
#include <boost/format.hpp>

#include "varoom/sam/bam_reader.hpp"
#include "varoom/sam/sam_flags.hpp"
#include "varoom/seq/gff.hpp"
#include "varoom/util/blob.hpp"
#include "varoom/util/files.hpp"
#include "varoom/util/ranges.hpp"

using namespace std;
using namespace boost;
using namespace varoom;
using namespace varoom::seq;

namespace std
{
    template<> struct hash<std::pair<size_t,size_t>>
    {
        std::size_t operator()(std::pair<size_t,size_t> const& x) const noexcept
        {
            std::size_t h1 = 0x3b7d530f5d27b987ULL + x.first * 0x317371dc620c3879ULL;
            std::size_t h2 = 0x1fce11d96fa69a03ULL + x.second * 0xfc3c906111de6791ULL;
            return h1 ^ h2;
        }
    };

    void print_row(std::ostream& p_out, const std::string& p_addr, size_t p_count, const std::string& p_ex, const map<string,vector<string>>& p_pmap)
    {
        vector<pair<string,string>> stk;
        stk.push_back(pair<string,string>(p_ex, p_ex));
        while (stk.size())
        {
            pair<string,string> itm = stk.back();
            stk.pop_back();
            auto itr = p_pmap.find(itm.first);
            if (itr == p_pmap.end())
            {
                p_out << p_addr << '\t' << p_count << '\t' << itm.second << std::endl;
                continue;
            }
            for (auto jtr = itr->second.begin(); jtr != itr->second.end(); ++jtr)
            {
                pair<string,string> jtm(*jtr, itm.second);
                jtm.second += "\t";
                jtm.second += *jtr;
                stk.push_back(jtm);
            }
        }
    }
}
// namespace std

namespace // anonymous
{
    using namespace std::chrono;
    using timer = std::chrono::high_resolution_clock;

    struct annotation
    {
        std::string chrom;
        uint64_t start1;
        uint64_t end1;
        std::string label;
    };

    class annotate_command : public varoom::command
    {
    public:
        annotate_command(const string& p_input_filename, const string& p_annot_filename,
                         const string& p_save_filename, const bool& p_load, const string& p_output_filename)
            : m_input_filename(p_input_filename),
              m_annot_filename(p_annot_filename),
              m_save_filename(p_save_filename),
              m_load(p_load),
              m_output_filename(p_output_filename)
        {
        }

        virtual void operator()()
        {
            if (m_load && m_save_filename.size() > 0)
            {
                throw std::runtime_error("cannot save and load in the one command invocation");
            }

            auto T0 = timer::now();
            map<string,map<ranges::range,vector<string>>> lmap;
            map<string,ranges> rmap;
            map<string,vector<string>> pmap;

            if (m_load)
            {
                input_file_holder_ptr inp = files::in(m_annot_filename);

                blob::load(**inp, [&](std::istream& p_in) mutable {
                    nlohmann::json jj;
                    p_in >> jj;
                    for (nlohmann::json::iterator itr = jj.begin(); itr != jj.end(); ++itr)
                    {
                        string k = itr.key();
                        nlohmann::json v = itr.value();
                        map<ranges::range,vector<string>> m;
                        for (nlohmann::json::iterator jtr = v.begin(); jtr != v.end(); ++jtr)
                        {
                            ranges::range r((*jtr)[0][0], (*jtr)[0][1]);
                            vector<string> s = (*jtr)[1];
                            m[r] = s;
                        }
                        lmap[k] = m;
                        rmap[k] = ranges();
                    }
                });

                blob::load(**inp, [&](std::istream& p_in) mutable {
                    nlohmann::json jj;
                    p_in >> jj;
                    for (nlohmann::json::iterator itr = jj.begin(); itr != jj.end(); ++itr)
                    {
                        string k = itr.key();
                        vector<string> v = itr.value();
                        pmap[k] = v;
                    }
                });
                for (auto itr = rmap.begin(); itr != rmap.end(); ++itr)
                {
                    itr->second.load(**inp);
                }
            }
            else
            {
                input_file_holder_ptr inp = files::in(m_annot_filename);
                const std::string kind("exon");
                map<std::string,ranges_builder> bldrs;

                bool found = false;
                annotation a;
                auto f = [&](const std::string& p_seqid, const std::string& p_source, const std::string& p_type,
                             const uint64_t& p_start1, const uint64_t& p_end1, const std::string& p_score,
                             const std::string& p_strand, const std::string& p_phase,
                             const std::vector<attr_and_val>& p_attributes) mutable {
                    found = false;

                    if (p_type != "exon" && p_type != "transcript" && p_type != "mRNA" && p_type != "gene")
                    {
                        return;
                    }

                    string id_str;
                    string name_str;
                    string parent_str;
                    for (size_t i = 0; i < p_attributes.size(); ++i)
                    {
                        if (p_attributes[i].first == "Name")
                        {
                            name_str = p_attributes[i].second;
                        }
                        if (p_attributes[i].first == "ID")
                        {
                            id_str = p_attributes[i].second;
                        }
                        if (p_attributes[i].first == "Parent")
                        {
                            parent_str = p_attributes[i].second;
                        }
                    }
                    if (id_str.size() == 0 && name_str.size() > 0)
                    {
                        id_str = name_str;
                    }
                    if (id_str.size() > 0)
                    {
                        if (parent_str.size() > 0)
                        {
                            pmap[id_str].push_back(parent_str);
                        }
                    }
                    if (id_str.size() > 0 && p_type == kind)
                    {
                        a.chrom = p_seqid;
                        a.start1 = p_start1;
                        a.end1 = p_end1;
                        a.label = id_str;
                        found = true;
                    }
                };

                size_t n = 0;
                size_t m = 0;
                gff3_reader r(**inp);
                while (r.next(f))
                {
                    ++n;
                    if (found == true)
                    {
                        ++m;
                        ranges_builder::range g(a.start1, a.end1);
                        bldrs[a.chrom].insert(g);
                        lmap[a.chrom][g].push_back(a.label);
                    }
                }
                std::cerr << "processed entries: " << n << std::endl;
                std::cerr << "matching entries: " << m << std::endl;
                std::cout << nlohmann::json(pmap) << std::endl;
                for (auto itr = bldrs.begin(); itr != bldrs.end(); ++itr)
                {
                    rmap[itr->first].make(itr->second);
                }

                if (m_save_filename.size() > 0)
                {
                    output_file_holder_ptr outp = files::out(m_save_filename);
                    blob::save(**outp, [&](std::ostream& p_out) {
                        nlohmann::json jj = nlohmann::json::object();
                        for (auto itr = lmap.begin(); itr != lmap.end(); ++itr)
                        {
                            nlohmann::json kk = nlohmann::json::array();
                            for (auto jtr = itr->second.begin(); jtr != itr->second.end(); ++jtr)
                            {
                                nlohmann::json pp({jtr->first.first, jtr->first.second});
                                kk.push_back(nlohmann::json({pp, jtr->second}));
                            }
                            jj[itr->first] = kk;
                        }
                        p_out << jj << std::endl;
                    });

                    blob::save(**outp, [&](std::ostream& p_out) {
                        nlohmann::json jj = pmap;
                        p_out << jj << std::endl;
                    });

                    for (auto itr = rmap.begin(); itr != rmap.end(); ++itr)
                    {
                        itr->second.save(**outp);
                    }
                    return;
                }
            }

            unordered_map<std::string,unordered_map<ranges::range,size_t>> A;
            vector<ranges::range> hitVec;
            auto T1 = timer::now();
            size_t n = 0;
            for (bam_reader r(m_input_filename); r.more(); ++r)
            {
                ++n;
                const sam_alignment& aln = *r;
                if (sam_flags::is_unmapped(aln.flags))
                {
                    continue;
                }
                if (aln.tlen < 0 || aln.tlen > 4096)
                {
                    continue;
                }
                ranges::range rx(aln.pos, aln.pos + aln.seq.size());
                auto itr = rmap.find(aln.chr);
                if (itr == rmap.end())
                {
                    continue;
                }
                itr->second.ranges_at(rx.first, rx.second, hitVec);
                if (hitVec.size() != 0)
                {
                    for (size_t i = 0; i < hitVec.size(); ++i)
                    {
                        A[aln.chr][hitVec[i]] += 1;
                    }
                    if (0)
                    {
                        nlohmann::json jj = nlohmann::json::array();
                        for (size_t i = 0; i < hitVec.size(); ++i)
                        {
                            jj.push_back(nlohmann::json({hitVec[i].first, hitVec[i].second}));
                        }
                        std::cout << aln.chr << '\t' << aln.pos << '\t' << jj << std::endl;
                    }
                }
            }
            std::cerr << "processed bam entries: " << n << std::endl;
            auto T2 = timer::now();
            for (auto itr = lmap.begin(); itr != lmap.end(); ++itr)
            {
                auto iitr = A.find(itr->first);
                for (auto jtr = itr->second.begin(); jtr != itr->second.end(); ++jtr)
                {
                    size_t c = 0;
                    if (iitr != A.end())
                    {
                        auto jjtr = iitr->second.find(jtr->first);
                        if (jjtr != iitr->second.end())
                        {
                            c = jjtr->second;
                        }
                    }
                    for (auto ktr = jtr->second.begin(); ktr != jtr->second.end(); ++ktr)
                    {
                        string addr = str(format("%s:%d-%d") % itr->first % jtr->first.first % jtr->first.second);
                        print_row(std::cout, addr, c, *ktr, pmap);
                    }
                }
            }
            auto T3 = timer::now();
            std::cerr << "creating range set: " << (1e-3*duration_cast<milliseconds>(T1 - T0).count()) << std::endl;
            std::cerr << "scanning bam: " << (1e-3*duration_cast<milliseconds>(T2 - T1).count()) << std::endl;
            std::cerr << "writing output: " << (1e-3*duration_cast<milliseconds>(T3 - T2).count()) << std::endl;
        }

    private:
        const string m_input_filename;
        const string m_annot_filename;
        const string m_save_filename;
        const bool m_load;
        const string m_output_filename;
    };

    class annotate_factory : public command_factory
    {
    public:
        annotate_factory() {}

        virtual json parse(const command_options& p_global_opts, const command_parameters& p_globals,
                           const std::vector<std::string>& p_args) const
        {
            namespace po = boost::program_options;

            command_options opts("compile a compact genome reference");
            opts.add_options()
                ("help,h", "produce help message")
                ("input,i", po::value<string>()->default_value("-"), "input filename, defaults to '-' for stdin")
                ("annot,a", po::value<string>(), "filename for gff3 file with features to annotate")
                ("save,S", po::value<string>(), "filename to save the indexed annotations to")
                ("load,L", "if specified, use the annotation name to load pre-indexed annotations")
                ("output,o", po::value<string>()->default_value("-"), "output filename, defaults to '-' for stdout")
                ;

            po::positional_options_description pos;

            po::variables_map vm;
            po::parsed_options parsed
                = po::command_line_parser(p_args).options(opts).positional(pos).run();
            po::store(parsed, vm);
            po::notify(vm);

            json params;
            if (p_globals.count("help") || vm.count("help"))
            {
                cout << p_global_opts << endl;
                cout << opts << endl;
                return params;
            }

            params["input"] = vm["input"].as<string>();
            params["annot"] = vm["annot"].as<string>();
            if (vm.count("save"))
            {
                params["save"] = vm["save"].as<string>();
            }
            else
            {
                params["save"] = string("");
            }
            params["load"] = (vm.count("load") ? true : false);
            params["output"] = vm["output"].as<string>();

            return params;
        }

        virtual command_ptr create(const json& p_params) const
        {
            if (!p_params.is_object())
            {
                return command_ptr();
            }
            string input_fn = p_params["input"];
            string annot_fn = p_params["annot"];
            string save_fn = p_params["save"];
            bool load = p_params["load"];
            string output_fn = p_params["output"];
            return command_ptr(new annotate_command(input_fn, annot_fn, save_fn, load, output_fn));
        }
    };
    
    bool reg = command_factory::add("annotate", command_factory_ptr(new annotate_factory));
}
// namespace anonymous

