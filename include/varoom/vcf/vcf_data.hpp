#ifndef VAROOM_VCF_VCF_DATA_HPP
#define VAROOM_VCF_VCF_DATA_HPP

#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace varoom
{
    namespace vcf
    {
        enum vcf_data_type { INT, FLT, STR, CHR, VEC, MAP };

        class vcf_data;
        typedef std::map<std::string,vcf_data> vcf_data_map;
        typedef std::vector<vcf_data> vcf_data_vector;

        class vcf_data
        {
        public:
            vcf_data(const std::int64_t& p_val)
            {
                ptrRep = new std::int64_t(p_val);
                intRep |= INT;
            }

            vcf_data(const double& p_val)
            {
                ptrRep = new double(p_val);
                intRep |= FLT;
            }

            vcf_data(const std::string& p_val)
            {
                ptrRep = new std::string(p_val);
                intRep |= STR;
            }

            vcf_data(const char p_val)
            {
                ptrRep = new char(p_val);
                intRep |= CHR;
            }

            vcf_data(const vcf_data_vector p_val)
            {
                ptrRep = new vcf_data_vector(p_val);
                intRep |= VEC;
            }

            vcf_data(const vcf_data_map p_val)
            {
                ptrRep = new vcf_data_map(p_val);
                intRep |= MAP;
            }

            vcf_data(const std::vector<vcf_data>& p_val)
            {
                ptrRep = new std::vector<vcf_data>(p_val);
                intRep |= VEC;
            }

            vcf_data_type type() const
            {
                return static_cast<vcf_data_type>(intRep & 0x7);
            }

            const std::int64_t& asInt() const
            {
                return *reinterpret_cast<const std::int64_t*>(detag());
            }

            const double& asFlt() const
            {
                return *reinterpret_cast<const double*>(detag());
            }

            const std::string& asStr() const
            {
                return *reinterpret_cast<const std::string*>(detag());
            }

            const char& asChr() const
            {
                return *reinterpret_cast<const char*>(detag());
            }

            const vcf_data_map& asMap() const
            {
                return *reinterpret_cast<const vcf_data_map*>(detag());
            }

            const vcf_data_vector& asVec() const
            {
                return *reinterpret_cast<const vcf_data_vector*>(detag());
            }

            // Copy constructor
            vcf_data(const vcf_data& p_other)
            {
                ptrRep = copy(p_other);
                intRep |= static_cast<std::intptr_t>(p_other.type());
            }

            // Assignment
            //
            vcf_data& operator=(const vcf_data& p_other)
            {
                if (this == &p_other)
                {
                    return *this;
                }
                release();
                ptrRep = copy(p_other);
                intRep |= static_cast<std::intptr_t>(p_other.type());
                return *this;
            }

            // Destructor
            //
            ~vcf_data()
            {
                release();
            }

        private:
            void release()
            {
                if (intRep == 0)
                {
                    return;
                }
                switch (type())
                {
                    case INT:
                        delete reinterpret_cast<const std::int64_t*>(detag());
                        break;
                    case FLT:
                        delete reinterpret_cast<const double*>(detag());
                        break;
                    case STR:
                        delete reinterpret_cast<const std::string*>(detag());
                        break;
                    case CHR:
                        delete reinterpret_cast<const char*>(detag());
                        break;
                    case MAP:
                        delete reinterpret_cast<const vcf_data_map*>(detag());
                        break;
                    case VEC:
                        delete reinterpret_cast<const vcf_data_vector*>(detag());
                        break;
                }
            }

            template <typename T>
            static void* repl(const T& p_x)
            {
                return new T(p_x);
            }

            static const void* copy(const vcf_data& p_other)
            {
                switch (p_other.type())
                {
                    case INT:
                        return repl(p_other.asInt());
                    case FLT:
                        return repl(p_other.asFlt());
                    case STR:
                        return repl(p_other.asStr());
                    case CHR:
                        return repl(p_other.asChr());
                    case MAP:
                        return repl(p_other.asMap());
                    case VEC:
                        return repl(p_other.asVec());
                }
                std::abort();
            }

            const void* detag() const
            {
                return reinterpret_cast<const void*>(intRep & ~static_cast<std::intptr_t>(0x7));
            }

            union {
                std::intptr_t intRep;
                const void*   ptrRep;
            };
        };
    }
    // namespace vcf
}
// namespace varoom

#endif // VAROOM_VCF_VCF_DATA_HPP
