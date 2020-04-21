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
        annotate_command(const string& p_input_filename, const string& p_annot_filename, const string& p_output_filename)
            : m_input_filename(p_input_filename),
              m_annot_filename(p_annot_filename),
              m_output_filename(p_output_filename)
        {
        }

        virtual void operator()()
        {
            input_file_holder_ptr inp = files::in(m_annot_filename);

            map<string,map<ranges::range,vector<string>>> lmap;
            map<string,ranges> rmap;
            {
                const std::string kind("gene");
                map<std::string,ranges_builder> bldrs;

                bool found = false;
                annotation a;
                auto f = [&](const std::string& p_seqid, const std::string& p_source, const std::string& p_type,
                             const uint64_t& p_start1, const uint64_t& p_end1, const std::string& p_score,
                             const std::string& p_strand, const std::string& p_phase,
                             const std::vector<attr_and_val>& p_attributes) mutable {
                    if (p_type == kind)
                    {
                        a.chrom = p_seqid;
                        a.start1 = p_start1;
                        a.end1 = p_end1;
                        a.label.clear();
                        for (size_t i = 0; i < p_attributes.size(); ++i)
                        {
                            if (p_attributes[i].first == "Name")
                            {
                                a.label = p_attributes[i].second;
                                break;
                            }
                        }
                        found = true;
                    }
                    else
                    {
                        found = false;
                    }
                };

                gff3_reader r(**inp);
                while (r.next(f))
                {
                    if (found == true)
                    {
                        ranges_builder::range g(a.start1, a.end1);
                        bldrs[a.chrom].insert(g);
                        lmap[a.chrom][g].push_back(a.label);
                    }
                }
                for (auto itr = bldrs.begin(); itr != bldrs.end(); ++itr)
                {
                    rmap[itr->first].make(itr->second);
                }
            }

            if (0)
            {
                output_file_holder_ptr outp = files::out(m_output_filename);
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
                for (auto itr = rmap.begin(); itr != rmap.end(); ++itr)
                {
                    itr->second.save(**outp);
                }
            }

            size_t hits = 0;
            double hitTime = 0;
            size_t misses = 0;
            double missTime = 0;

            unordered_map<std::string,unordered_map<ranges::range,size_t>> A;
            vector<ranges::range> hitVec;
            auto T0 = timer::now();
            for (bam_reader r(m_input_filename); r.more(); ++r)
            {
                const sam_alignment& aln = *r;
                if (sam_flags::is_unmapped(aln.flags))
                {
                    continue;
                }
                if (aln.tlen < 0 || aln.tlen > 4096)
                {
                    continue;
                }
                ranges::range rx(aln.pos, aln.pos + aln.tlen);
                auto itr = rmap.find(aln.chr);
                if (itr == rmap.end())
                {
                    continue;
                }
                auto t0 = timer::now();
                itr->second.ranges_at(rx.first, rx.second, hitVec);
                auto t1 = timer::now();
                std::chrono::duration<double> d = t1 - t0;
                if (hitVec.size() != 0)
                {
                    hits += 1;
                    hitTime += d.count();
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
                else
                {
                    misses += 1;
                    missTime += d.count();
                }
            }
            auto T1 = timer::now();
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
                        std::cout << str(format("%s\t%s:%d-%d\t%d") % (*ktr) % itr->first % jtr->first.first % jtr->first.second % c) << std::endl;
                    }
                }
            }
        }

    private:
        const string m_input_filename;
        const string m_annot_filename;
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
            string output_fn = p_params["output"];
            return command_ptr(new annotate_command(input_fn, annot_fn, output_fn));
        }
    };
    
    bool reg = command_factory::add("annotate", command_factory_ptr(new annotate_factory));
}
// namespace anonymous

