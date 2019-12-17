#ifndef VAROOM_UTIL_GZIP_STRINGS_HPP
#define VAROOM_UTIL_GZIP_STRINGS_HPP

#include <sstream>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

namespace varoom
{
    class gzip_strings {
        public:
            static void compress(const std::string& p_src, std::string& p_dst)
            {
                std::stringstream compressed;
                std::stringstream origin(p_src);

                boost::iostreams::filtering_streambuf<boost::iostreams::input> out;
                out.push(boost::iostreams::gzip_compressor(boost::iostreams::gzip_params(boost::iostreams::gzip::best_compression)));
                out.push(origin);
                boost::iostreams::copy(out, compressed);

                p_dst = compressed.str();
            }

            static std::string compress(const std::string& p_src)
            {
                std::string res;
                compress(p_src, res);
                return res;
            }

            static void decompress(const std::string& p_src, std::string& p_dst)
            {
                std::stringstream compressed(p_src);
                std::stringstream decompressed;

                boost::iostreams::filtering_streambuf<boost::iostreams::input> out;
                out.push(boost::iostreams::gzip_decompressor());
                out.push(compressed);
                boost::iostreams::copy(out, decompressed);

                p_dst = decompressed.str();
            }

            static std::string decompress(const std::string& p_src)
            {
                std::string res;
                decompress(p_src, res);
                return res;
            }
    };
}
// namespace varoom

#endif // VAROOM_UTIL_GZIP_STRINGS_HPP
