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
            static std::string compress(const std::string& p_src)
            {
                std::stringstream compressed;
                std::stringstream origin(p_src);

                boost::iostreams::filtering_streambuf<boost::iostreams::input> out;
                out.push(boost::iostreams::gzip_compressor(boost::iostreams::gzip_params(boost::iostreams::gzip::best_compression)));
                out.push(origin);
                boost::iostreams::copy(out, compressed);

                return compressed.str();
            }

            static std::string decompress(const std::string& p_src)
            {
                std::stringstream compressed(p_src);
                std::stringstream decompressed;

                boost::iostreams::filtering_streambuf<boost::iostreams::input> out;
                out.push(boost::iostreams::gzip_decompressor());
                out.push(compressed);
                boost::iostreams::copy(out, decompressed);

                return decompressed.str();
            }
    };
}
// namespace varoom

#endif // VAROOM_UTIL_GZIP_STRINGS_HPP
