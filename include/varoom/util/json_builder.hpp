#ifndef VAROOM_UTIL_JSON_BUILDER_HPP
#define VAROOM_UTIL_JSON_BUILDER_HPP

#include <nlohmann/json.hpp>

namespace varoom
{
    class json_writer
    {
    private:
        enum state { OK, DONE, ARRAY_START, IN_ARRAY, OBJECT_START, IN_OBJECT, SEEN_KEY };

    public:
        json_writer(std::ostream& p_out)
            : m_out(p_out)
        {
            m_stack.push_back(OK);
        }

        void null()
        {
            write_value([](std::ostream& p_out) {
                p_out << "null";
            });
        }

        void boolean(bool p_bool)
        {
            write_value([=](std::ostream& p_out) {
                p_out << (p_bool ? "true" : "false");
            });
        }

        void number_integer(int64_t p_num)
        {
            write_value([=](std::ostream& p_out) {
                p_out << p_num;
            });
        }

        void number_unsigned(uint64_t p_num)
        {
            write_value([=](std::ostream& p_out) {
                p_out << p_num;
            });
        }

        void number_float(double p_num)
        {
            write_value([=](std::ostream& p_out) {
                p_out << p_num;
            });
        }

        void string(const std::string& p_str)
        {
            write_value([=](std::ostream& p_out) {
                p_out << nlohmann::json(p_str);
            });
        }

        void begin_array()
        {
            switch (m_stack.back())
            {
                case OK:
                {
                    m_out << '[';
                    m_stack.push_back(ARRAY_START);
                    break;
                }
                case IN_ARRAY:
                {
                    m_out << ", ";
                }
                case ARRAY_START:
                {
                    m_out << '[';
                    m_stack.push_back(ARRAY_START);
                    break;
                };
                case SEEN_KEY:
                {
                    m_out << ":[";
                    m_stack.push_back(ARRAY_START);
                    break;
                }
                default:
                {
                    throw std::runtime_error("cannot write array in this context.");
                }
            }
        }

        void end_array()
        {
            switch (m_stack.back())
            {
                case ARRAY_START:
                case IN_ARRAY:
                {
                    m_out << ']';
                    break;
                }
                default:
                {
                    throw std::runtime_error("cannot end array in this context.");
                }
            }
            m_stack.pop_back();
            switch (m_stack.back())
            {
                case OK:
                {
                    m_stack.back() = DONE;
                    break;
                }
                case ARRAY_START:
                case IN_ARRAY:
                {
                    m_stack.back() = IN_ARRAY;
                    break;
                }
                case SEEN_KEY:
                {
                    m_stack.back() = IN_OBJECT;
                    break;
                }
                default:
                {
                    throw std::runtime_error("internal error after array.");
                }
            }
        }

        void begin_object()
        {
            switch (m_stack.back())
            {
                case OK:
                {
                    m_out << '{';
                    m_stack.push_back(OBJECT_START);
                    break;
                }
                case IN_ARRAY:
                {
                    m_out << ", ";
                }
                case ARRAY_START:
                {
                    m_out << '{';
                    m_stack.push_back(OBJECT_START);
                    break;
                };
                case SEEN_KEY:
                {
                    m_out << ":{";
                    m_stack.push_back(OBJECT_START);
                    break;
                }
                default:
                {
                    throw std::runtime_error("cannot write object in this context.");
                }
            }
        }

        void end_object()
        {
            switch (m_stack.back())
            {
                case OBJECT_START:
                case IN_OBJECT:
                {
                    m_out << '}';
                    break;
                }
                default:
                {
                    throw std::runtime_error("cannot end object in this context.");
                }
            }
            m_stack.pop_back();
            switch (m_stack.back())
            {
                case OK:
                {
                    m_stack.back() = DONE;
                    break;
                }
                case ARRAY_START:
                case IN_ARRAY:
                {
                    m_stack.back() = IN_ARRAY;
                    break;
                }
                case SEEN_KEY:
                {
                    m_stack.back() = IN_OBJECT;
                    break;
                }
                default:
                {
                    throw std::runtime_error("internal error after object.");
                }
            }
        }

        void key(const std::string& p_key)
        {
            switch (m_stack.back())
            {
                case OBJECT_START:
                case IN_OBJECT:
                {
                    m_out << nlohmann::json(p_key);
                    m_stack.back() = SEEN_KEY;
                    break;
                }
                default:
                {
                    throw std::runtime_error("cannot end object in this context.");
                }
            }
        }

    private:
        template <class F>
        void write_value(F p_func)
        {
            switch (m_stack.back())
            {
                case OK:
                {
                    p_func(m_out);
                    m_stack.back() = DONE;
                    break;
                }
                case IN_ARRAY:
                {
                    m_out << ", ";
                }
                case ARRAY_START:
                {
                    p_func(m_out);
                    m_stack.back() = IN_ARRAY;
                    break;
                }
                case SEEN_KEY:
                {
                    m_out << ':';
                    p_func(m_out);
                    m_stack.back() = IN_OBJECT;
                    break;
                }
                default:
                {
                    throw std::runtime_error("cannot place NULL in this context.");
                }
            }
        }

        std::ostream& m_out;
        std::vector<state> m_stack;
    };
}
// namespace varoom

#endif // VAROOM_UTIL_JSON_BUILDER_HPP
