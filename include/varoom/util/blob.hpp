#ifndef VAROOM_UTIL_BLOB_HPP
#define VAROOM_UTIL_BLOB_HPP

#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/stream.hpp>

namespace varoom
{
    class blob
    {
    public:
        static void save(std::ostream& p_out, std::function<void(std::ostream&)> p_saver)
        {
            std::vector<char> b;
            boost::iostreams::back_insert_device<std::vector<char>> d(b);
            boost::iostreams::stream<boost::iostreams::back_insert_device<std::vector<char>>> s(d);
            p_saver(s);
            s.flush();

            uint64_t z = b.size();
            p_out.write(reinterpret_cast<const char*>(&z), sizeof(uint64_t));
            p_out.write(b.data(), z);
        }

        static void save_gz(std::ostream& p_out, std::function<void(std::ostream&)> p_saver)
        {
            std::vector<char> b;
            boost::iostreams::back_insert_device<std::vector<char>> dev(b);
            boost::iostreams::filtering_ostream out;
            out.push(boost::iostreams::gzip_compressor(boost::iostreams::gzip_params(boost::iostreams::gzip::best_compression)));
            out.push(dev);
            p_saver(out);
            out.flush();

            uint64_t z = b.size();
            p_out.write(reinterpret_cast<const char*>(&z), sizeof(uint64_t));
            p_out.write(b.data(), z);
        }

        static void load(std::istream& p_in, std::function<void(std::istream&)> p_loader)
        {
            uint64_t z;
            p_in.read(reinterpret_cast<char*>(&z), sizeof(uint64_t));
            std::string s;
            s.resize(z);
            p_in.read(&s[0], z);
            boost::iostreams::array_source dd(&s[0], z);
            boost::iostreams::stream<boost::iostreams::array_source> ss(dd);
            p_loader(ss);
        }

        static void load_gz(std::istream& p_in, std::function<void(std::istream&)> p_loader)
        {
            uint64_t z;
            p_in.read(reinterpret_cast<char*>(&z), sizeof(uint64_t));
            std::string s;
            s.resize(z);
            p_in.read(&s[0], z);

            boost::iostreams::array_source dev(&s[0], z);
            boost::iostreams::filtering_istream in;
            in.push(boost::iostreams::gzip_decompressor());
            in.push(dev);
            p_loader(in);
        }
    };
}
// namespace varoom

#endif // VAROOM_UTIL_BLOB_HPP
