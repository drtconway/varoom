#include "varoom/command.hpp"
#include "varoom/util/files.hpp"
#include "varoom/util/subtext.hpp"
#include "varoom/util/table.hpp"
#include "varoom/klbam/types.hpp"

#include <boost/lexical_cast.hpp>
#include <boost/range/iterator_range.hpp>

using namespace std;
using namespace boost;
using namespace varoom;
using namespace varoom::klbam;

namespace // anonymous
{
    typedef vector<string> strings;

    template <int I>
    void add(const counts::tuple& p_lhs, const counts::tuple& p_rhs, counts::tuple& p_res)
    {
        std::get<I>(p_res) = std::get<I>(p_lhs) + std::get<I>(p_rhs);
    }

    void parse_other(const string& p_str, map<string,size_t>& p_res)
    {
        subtext s(p_str);
        vector<subtext> ps;
        vector<subtext> qs;
        s.split(';', ps);
        for (size_t i = 0; i < ps.size(); ++i)
        {
            ps[i].split('=', qs);
            if (qs.size() != 2)
            {
                throw std::runtime_error("bad other string: " + p_str);
            }
            p_res[qs[0]] += boost::lexical_cast<size_t>(boost::make_iterator_range(qs[1].first, qs[1].second));
        }
    }

    string serialize_other(const map<string,size_t>& p_oth)
    {
        string res;
        for (auto itr = p_oth.begin(); itr != p_oth.end(); ++itr)
        {
            if (res.size() > 0)
            {
                res.push_back(';');
            }
            string s = itr->first + string("=") + boost::lexical_cast<string>(itr->second);
            res.insert(res.end(), s.begin(), s.end());
        }
        return res;
    }

    string merge_other(const string& p_lhs, const string& p_rhs)
    {
        if (p_rhs.size() == 0)
        {
            return p_lhs;
        }
        if (p_lhs.size() == 0)
        {
            return p_rhs;
        }

        map<string,size_t> m;
        parse_other(p_lhs, m);
        parse_other(p_rhs, m);
        return serialize_other(m);
    }

    bool row_less(const counts::tuple& p_lhs, const counts::tuple& p_rhs)
    {
        if (std::get<counts::chr>(p_lhs) < std::get<counts::chr>(p_rhs))
        {
            return true;
        }
        if (std::get<counts::chr>(p_rhs) < std::get<counts::chr>(p_lhs))
        {
            return false;
        }
        return std::get<counts::pos>(p_lhs) < std::get<counts::pos>(p_rhs);
    }

    void merge_rows(const counts::tuple& p_lhs, const counts::tuple& p_rhs, counts::tuple& p_res)
    {
        std::get<counts::chr>(p_res) = std::get<counts::chr>(p_lhs);
        std::get<counts::pos>(p_res) = std::get<counts::pos>(p_lhs);
        add<counts::nA>(p_lhs, p_rhs, p_res);
        add<counts::nC>(p_lhs, p_rhs, p_res);
        add<counts::nG>(p_lhs, p_rhs, p_res);
        add<counts::nT>(p_lhs, p_rhs, p_res);
        add<counts::nN>(p_lhs, p_rhs, p_res);
        add<counts::nX0>(p_lhs, p_rhs, p_res);
        add<counts::nX1>(p_lhs, p_rhs, p_res);
        add<counts::nI>(p_lhs, p_rhs, p_res);
        add<counts::nD>(p_lhs, p_rhs, p_res);
        std::get<counts::oth>(p_res) = merge_other(std::get<counts::oth>(p_lhs), std::get<counts::oth>(p_rhs));
    }

    class merge_counts_command : public varoom::command
    {
    public:
        merge_counts_command(const strings& p_input_filenames,
                            const string& p_output_filename)
            : m_input_filenames(p_input_filenames),
              m_output_filename(p_output_filename)
        {
        }

        virtual void operator()()
        {
            std::function<bool(const counts::tuple&,const counts::tuple&)> l = row_less;
            std::function<void(const counts::tuple&,const counts::tuple&,counts::tuple&)> f = merge_rows;

            counts::table tbl;
            counts::table tmp;

            for (size_t n = 0; n < m_input_filenames.size(); ++n)
            {
                tmp.clear();
                input_file_holder_ptr inp = files::in(m_input_filenames[n]);
                counts::istream_reader lhs(**inp, true);
                counts::read_iterator rhs(tbl);
                counts::write_iterator res(tmp);
                table_utils::merge(lhs, rhs, res, l, f);
                std::swap(tmp, tbl);
            }

            output_file_holder_ptr outp = files::out(m_output_filename);
            counts::write(**outp, tbl, counts::labels());
        }

    private:
        const strings m_input_filenames;
        const string m_output_filename;
    };

    class merge_counts_factory : public command_factory
    {
    public:
        merge_counts_factory() {}

        virtual json parse(const command_options& p_global_opts, const command_parameters& p_globals,
                           const std::vector<std::string>& p_args) const
        {
            namespace po = boost::program_options;

            command_options opts("merge counts");
            opts.add_options()
                ("help,h", "produce help message")
                ("inputs", po::value<strings>(), "input filename, defaults to '-' for stdin")
                ("output,o", po::value<string>()->default_value("-"), "output filename, defaults to '-' for stdout")
                ;

            po::positional_options_description pos;
            pos.add("inputs", -1);

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

            if (vm.count("inputs"))
            {
                params["inputs"] = vm["inputs"].as<strings>();
            }
            else
            {
                strings ss{"-"};
                params["inputs"] = ss;
            }
            params["output"] = vm["output"].as<string>();

            return params;
        }

        virtual command_ptr create(const json& p_params) const
        {
            if (!p_params.is_object())
            {
                return command_ptr();
            }
            strings input_fns = p_params["inputs"];
            string output_fn = p_params["output"];
            return command_ptr(new merge_counts_command(input_fns, output_fn));
        }
    };

    bool reg = command_factory::add("merge", command_factory_ptr(new merge_counts_factory));
}
// namespace anonymous
