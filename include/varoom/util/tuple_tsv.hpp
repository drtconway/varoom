#ifndef VAROOM_UTIL_TUPLE_TSV_HPP
#define VAROOM_UTIL_TUPLE_TSV_HPP

#ifndef VAROOM_UTIL_SUBTEXT_HPP
#include "varoom/util/subtext.hpp"
#endif

#include <boost/any.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/range/iterator_range.hpp>

namespace varoom
{
    namespace detail
    {
        template<std::size_t I = 0, typename... Tp>
        typename std::enable_if<I == sizeof...(Tp), void>::type
        to_tuple(std::vector<varoom::subtext>&, std::tuple<Tp...> &)
        {
        }

        template<std::size_t I = 0, typename... Tp>
        typename std::enable_if<I < sizeof...(Tp), void>::type
        to_tuple(std::vector<varoom::subtext>& ss, std::tuple<Tp...> & t)
        {
            varoom::subtext s = ss[I];
            auto r = boost::make_iterator_range(s.first, s.second);
            std::get<I>(t) = boost::lexical_cast<typename std::tuple_element<I, std::tuple<Tp...>>::type>(r);
            to_tuple<I + 1, Tp...>(ss, t);
        }

        template<std::size_t I = 0, typename... Tp>
        typename std::enable_if<I == sizeof...(Tp), void>::type
        to_tuple(std::vector<varoom::subtext>&, std::tuple<Tp...> &, const std::vector<boost::any>&)
        {
        }

        template<std::size_t I = 0, typename... Tp>
        typename std::enable_if<I < sizeof...(Tp), void>::type
        to_tuple(std::vector<varoom::subtext>& ss, std::tuple<Tp...> & t, const std::vector<boost::any>& cs)
        {
            typedef typename std::tuple_element<I, std::tuple<Tp...>>::type elem_type;
            typedef std::function<elem_type(const varoom::subtext&)> func_type;
            varoom::subtext s = ss[I];
            func_type f = boost::any_cast<func_type>(cs[I]);
            std::get<I>(t) = f(s);
            to_tuple<I + 1, Tp...>(ss, t, cs);
        }

    }
    // namespace detail

    class tuple_tsv
    {
    public:
        template<typename... Tp>
        static void to_tuple(std::vector<varoom::subtext>& ss, std::tuple<Tp...> & t)
        {
            detail::to_tuple(ss, t);
        }

        template<typename... Tp>
        static void to_tuple(std::vector<varoom::subtext>& ss, std::tuple<Tp...> & t, const std::vector<boost::any>& cs)
        {
            detail::to_tuple(ss, t, cs);
        }

        static std::function<std::string(const varoom::subtext&)> to_str()
        {
            return to_str_impl;
        }

        static std::function<int(const varoom::subtext&)> to_int()
        {
            return to_X_impl<int>;
        }

        static std::function<double(const varoom::subtext&)> to_flt()
        {
            return to_X_impl<double>;
        }

        static std::function<uint32_t(const varoom::subtext&)> to_uint32()
        {
            return to_X_impl<uint32_t>;
        }

        static std::function<int32_t(const varoom::subtext&)> to_int32()
        {
            return to_X_impl<int32_t>;
        }

        static std::function<uint64_t(const varoom::subtext&)> to_uint64()
        {
            return to_X_impl<uint64_t>;
        }

        static std::function<int64_t(const varoom::subtext&)> to_int64()
        {
            return to_X_impl<int64_t>;
        }

        static std::function<size_t(const varoom::subtext&)> to_size_t()
        {
            return to_X_impl<size_t>;
        }

    private:
        static std::string to_str_impl(const varoom::subtext& p_st)
        {
            return p_st;
        }

        template <typename T>
        static T to_X_impl(const varoom::subtext& p_st)
        {
            auto r = boost::make_iterator_range(p_st.first, p_st.second);
            return boost::lexical_cast<T>(r);
        }
    };
}
// namespace varoom

#endif // VAROOM_UTIL_TUPLE_TSV_HPP
