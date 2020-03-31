#include "varoom/command.hpp"

#include "varoom/seq/fastq.hpp"
#include "varoom/util/files.hpp"

#include <condition_variable>
#include <mutex>
#include <queue>
#include <thread>
#include <boost/format.hpp>

using namespace std;
using namespace boost;
using namespace varoom;
using namespace varoom::seq;

namespace // anonymous
{
    using namespace std::chrono;
    using strings = vector<string>;
    using timer = std::chrono::high_resolution_clock;

    map<size_t,size_t> counts;

    class read_multibox
    {
    public:
        static constexpr size_t max_size = 1000;

        read_multibox()
            : m_done(false)
        {
        }

        ~read_multibox()
        {
            if (0)
            {
                for (auto itr = m_counts.begin(); itr != m_counts.end(); ++itr)
                {
                    counts[itr->first] += itr->second;
                }
            }
        }

        void put(const fastq_read& p_read)
        {
            std::unique_lock<std::mutex> lk(m_mut);
            while (m_reads.size() >= max_size)
            {
                m_cond.wait(lk);
            }
            m_reads.push_back(p_read);
            lk.unlock();
            m_cond.notify_one();
        }

        void done()
        {
            std::unique_lock<std::mutex> lk(m_mut);
            while (m_reads.size() > 0)
            {
                m_cond.wait(lk);
            }
            m_done = true;
            lk.unlock();
            m_cond.notify_one();
        }

        bool get(vector<fastq_read>& p_reads)
        {
            std::unique_lock<std::mutex> lk(m_mut);
            while (!m_done && m_reads.size() == 0)
            {
                m_cond.wait(lk);
            }
            if (m_reads.size() > 0)
            {
                if (0)
                {
                    m_counts[m_reads.size()] += 1;
                }
                m_reads.swap(p_reads);
                m_reads.clear();
                lk.unlock();
                m_cond.notify_one();
                return true;
            }
            return false;
        }

    private:
        std::mutex m_mut;
        std::condition_variable m_cond;
        bool m_done;
        std::vector<fastq_read> m_reads;
        std::map<size_t,size_t> m_counts;
    };

    class read_box
    {
    public:
        read_box()
            : m_got_one(false), m_done(false)
        {
        }

        void put(const fastq_read& p_read)
        {
            std::unique_lock<std::mutex> lk(m_mut);
            while (m_got_one)
            {
                m_cond.wait(lk);
            }
            m_got_one = true;
            m_read = p_read;
            lk.unlock();
            m_cond.notify_one();
        }

        void done()
        {
            std::unique_lock<std::mutex> lk(m_mut);
            while (m_got_one)
            {
                m_cond.wait(lk);
            }
            m_done = true;
            lk.unlock();
            m_cond.notify_one();
        }

        bool get(fastq_read& p_read)
        {
            std::unique_lock<std::mutex> lk(m_mut);
            while (!m_done && !m_got_one)
            {
                m_cond.wait(lk);
            }
            if (m_got_one)
            {
                p_read = m_read;
                m_got_one = false;
                lk.unlock();
                m_cond.notify_one();
                return true;
            }
            return false;
        }

    private:
        std::mutex m_mut;
        std::condition_variable m_cond;
        bool m_got_one;
        bool m_done;
        fastq_read m_read;
    };

    class background_fastq_writer
    {
    public:
        background_fastq_writer(const std::string& p_filename)
            : m_output_holder(files::out(p_filename)), m_output(**m_output_holder)
        {
 #if 0
            m_thread = std::thread([this]() mutable {
                fastq_read r;
                while (m_box.get(r))
                {
                    fastq_writer::write(m_output, r);
                }
            });
#endif
            m_thread = std::thread([this]() mutable {
                vector<fastq_read> rs;
                while (m_box.get(rs))
                {
                    for (size_t i = 0; i < rs.size(); ++i)
                    {
                        fastq_writer::write(m_output, rs[i]);
                    }
                }
            });
        }

        ~background_fastq_writer()
        {
            m_box.done();
            m_thread.join();
        }

        void write(const fastq_read& p_read)
        {
            m_box.put(p_read);
        }

        void close()
        {
            m_box.done();
        }

    private:
        varoom::output_file_holder_ptr m_output_holder;
        std::ostream& m_output;
        read_multibox m_box;
        std::thread m_thread;
    };
    using background_fastq_writer_ptr = std::shared_ptr<background_fastq_writer>;

    class scatter_command : public varoom::command
    {
    public:
        scatter_command(const string& p_input_filename, const strings& p_output_filenames)
            : m_input_filename(p_input_filename),
              m_output_filenames(p_output_filenames)
        {
        }

#if 1
        virtual void operator()()
        {
            {
                vector<background_fastq_writer_ptr> outs;
                for (size_t i = 0; i < m_output_filenames.size(); ++i)
                {
                    outs.push_back(background_fastq_writer_ptr(new background_fastq_writer(m_output_filenames[i])));
                }
                const size_t N = outs.size();

                input_file_holder_ptr inp = files::in(m_input_filename);

                auto t0 = timer::now();
                size_t n = 0;
                for (fastq_reader r(**inp); r.more(); ++r, ++n)
                {
                    while (n >= N)
                    {
                        n -= N;
                    }
                    outs[n]->write(*r);
                }
                auto t1 = timer::now();
                double d = duration_cast<milliseconds>(t1-t0).count() / 1000.0;
                std::cerr << "elapsed time: " << d << std::endl;
            }
            if (0)
            {
                std::cerr << "block sizes:" << std::endl;
                for (auto itr = counts.begin(); itr != counts.end(); ++itr)
                {
                    std::cerr << itr->first << '\t' << itr->second << std::endl;
                }
            }
        }
#endif
#if 0
        virtual void operator()()
        {
            vector<output_file_holder_ptr> outs;
            for (size_t i = 0; i < m_output_filenames.size(); ++i)
            {
                outs.push_back(files::out(m_output_filenames[i]));
            }
            const size_t N = outs.size();

            input_file_holder_ptr inp = files::in(m_input_filename);

            auto t0 = timer::now();
            size_t n = 0;
            for (fastq_reader r(**inp); r.more(); ++r, ++n)
            {
                while (n >= N)
                {
                    n -= N;
                }
                fastq_writer::write(**(outs[n]), *r);
            }
            auto t1 = timer::now();
            double d = duration_cast<milliseconds>(t1-t0).count() / 1000.0;
            std::cerr << "elapsed time: " << d << std::endl;
        }
#endif

    private:
        const string m_input_filename;
        const strings m_output_filenames;
    };

    class scatter_factory : public command_factory
    {
    public:
        scatter_factory() {}

        virtual json parse(const command_options& p_global_opts, const command_parameters& p_globals,
                           const std::vector<std::string>& p_args) const
        {
            namespace po = boost::program_options;

            command_options opts("compile a compact genome reference");
            opts.add_options()
                ("help,h", "produce help message")
                ("input,i", po::value<string>()->default_value("-"), "input filename, defaults to '-' for stdin")
                ("directory,d", po::value<string>(), "directory for output files")
                ("number,n", po::value<int>()->default_value(2), "number of files to split the output into")
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

            int n = vm["number"].as<int>();
            string d = vm["directory"].as<string>();

            strings outputs;
            for (int i = 0; i < n; ++i)
            {
                outputs.push_back(str(format("%s/part-%03d.fastq.gz") % d % i));
            }

            params["input"] = vm["input"].as<string>();
            //params["output"] = vm["output"].as<strings>();
            params["output"] = outputs;

            return params;
        }

        virtual command_ptr create(const json& p_params) const
        {
            if (!p_params.is_object())
            {
                return command_ptr();
            }
            string input_fn = p_params["input"];
            strings output_fns = p_params["output"];
            return command_ptr(new scatter_command(input_fn, output_fns));
        }
    };
    
    bool reg = command_factory::add("scatter", command_factory_ptr(new scatter_factory));
}
// namespace anonymous

