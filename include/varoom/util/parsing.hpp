#ifndef VAROOM_UTIL_PARSING_HPP
#define VAROOM_UTIL_PARSING_HPP

#include <functional>
#include <memory>
#include <type_traits>

namespace varoom
{
    namespace parsing
    {
        struct unit {};

        namespace detail
        {
            template <typename T>
            class suspension
            {
            public:
                explicit suspension(std::function<std::unique_ptr<T>()> p_func)
                    : m_func(p_func)
                {
                }

                const T& get() const
                {
                    if (!m_memo)
                    {
                        m_memo = m_func();
                    }
                    return *m_memo;
                }

            private:
                std::function<std::unique_ptr<T>()> m_func;
                mutable std::unique_ptr<T> m_memo;
            };

            template <typename T>
            class stream;

            template <typename T>
            class cell
            {
            public:
                cell() {}

                cell(T p_head)
                    : m_head(p_head)
                {
                }

                cell(T p_head, stream<T> p_tail)
                    : m_head(p_head), m_tail(std::move(p_tail))
                {
                }

                const T& head() const
                {
                    return m_head;
                }

                stream<T> tail() const
                {
                    return m_tail;
                }

            private:
                T m_head;
                stream<T> m_tail;
            };

            template <typename T>
            class stream
            {
            public:
                stream() {}

                stream(std::function<std::unique_ptr<cell<T>>()> p_func)
                    : m_lazy_cell(std::make_shared<suspension<cell<T>>>(p_func))
                {
                }

                //stream(stream&& p_other)
                //    : m_lazy_cell(std::move(p_other.m_lazy_cell))
                //{
                //}

                //stream& operator=(stream&& p_other)
                //{
                //    m_lazy_cell = std::move(p_other.m_lazy_cell);
                //    return *this;
                //}

                bool empty() const
                {
                    return !m_lazy_cell;
                }

                const T& head() const
                {
                    return m_lazy_cell->get().head();
                }

                stream<T> tail() const
                {
                    return m_lazy_cell->get().tail();
                }

                stream take(size_t p_n) const
                {
                    if (p_n == 0 || empty())
                    {
                        return stream();
                    }
                    auto h = head();
                    auto t = tail();
                    return stream([=]() {
                        return std::unique_ptr<cell<T>>(new cell<T>(h, t.take(p_n - 1)));
                    });
                }

            private:
                std::shared_ptr<suspension<cell<T>>> m_lazy_cell;
            };

            template <typename T, typename F>
            void for_each(stream<T> p_stream, F p_func)
            {
                while (!p_stream.empty())
                {
                    p_func(p_stream.head());
                    p_stream = p_stream.tail();
                }
            }

            template <typename T, typename F>
            auto fmap(stream<T> p_stream, F p_func) -> stream<decltype(p_func(p_stream.head()))>
            {
                using U = decltype(p_func(p_stream.head()));
                static_assert(std::is_convertible<F, std::function<U(T)>>::value, "fmap requires a function type U(T)");

                if (p_stream.empty())
                {
                    return stream<U>();
                }

                return stream<U>([p_stream, p_func]() {
                    return std::unique_ptr<cell<U>>(new cell<U>(p_func(p_stream.head()), fmap(p_stream.tail(), p_func)));
                });
            }

            template <typename T, typename F>
            auto fmapv(stream<T> p_stream, F p_func) -> stream<decltype(p_func())>
            {
                using U = decltype(p_func());
                static_assert(std::is_convertible<F, std::function<U(T)>>::value, "fmap requires a function type U()");

                if (p_stream.empty())
                {
                    return stream<U>();
                }

                return stream<U>([p_stream, p_func]() {
                    return std::unique_ptr<cell<U>>(new cell<U>(p_func(), fmapv(p_stream.tail(), p_func)));
                });
            }

            template <typename T>
            stream<T> concat(stream<T> p_lhs, stream<T> p_rhs)
            {
                if (p_lhs.empty())
                {
                    return p_rhs;
                }
                return stream<T>([p_lhs, p_rhs]() {
                    return std::unique_ptr<cell<T>>(new cell<T>(p_lhs.head(), concat(p_lhs.tail(), p_rhs)));
                });
            }

            template <typename T>
            stream<T> mjoin(stream<stream<T>> p_streams)
            {
                while (!p_streams.empty() && p_streams.head().empty())
                {
                    p_streams = p_streams.tail();
                }

                if (p_streams.empty())
                {
                    return stream<T>();
                }

                return stream<T>([p_streams]() {
                    stream<T> h = p_streams.head();
                    return std::unique_ptr<cell<T>>(new cell<T>(h.head(), concat(h.tail(), mjoin(p_streams.tail()))));
                });
            }

            template <typename T, typename F>
            auto mbind(stream<T> p_stream, F p_func) -> decltype(p_func(p_stream.head()))
            {
                return mjoin(fmap(p_stream, p_func));
            }

            template <typename T, typename F>
            auto mthen(stream<T> p_stream, F p_func) -> decltype(p_func())
            {
                return mjoin(fmapv(p_stream, p_func));
            }

            template <typename T>
            stream<T> mreturn(T p_value)
            {
                return stream<T>([p_value]() {
                    return std::unique_ptr<cell<T>>(new cell<T>(p_value));
                });
            }

            inline stream<unit> guard(bool p_value)
            {
                if (p_value)
                {
                    return mreturn(unit());
                }
                else
                {
                    return stream<unit>();
                }
            }
        }
        // namespace detail

        template <typename Input, typename T>
        using parse_result = detail::stream<std::pair<T,Input>>;

        template <typename Input, typename T>
        using parser = std::function<parse_result<Input,T>(Input)>;
    }
    // namespace parsing
}
// namespace varoom

#endif // VAROOM_UTIL_PARSING_HPP
