#ifndef VAROOM_FUNDS_PARSING_HPP
#define VAROOM_FUNDS_PARSING_HPP

#ifndef VAROOM_FUNDS_LIST_HPP
#include "varoom/funds/list.hpp"
#endif

namespace varoom
{
    namespace funds
    {
        using state = std::string;

        template <typename S, typename T>
        using pair_list = varoom::funds::list<std::pair<T,S>>;

        template <typename S, typename T>
        using parser = std::function<pair_list<S,T>(S)>;

        template <typename S, typename T>
        pair_list<S,T> parse(parser<S,T> p, S s)
        {
            return p(s);
        }
        
        template <typename T, typename S>
        parser<S,T> yield(T t)
        {
            return [t] (S s) { return pair_list<S,T>(std::make_pair(t, s)); };
        }

        template <typename U, typename T, typename S, typename X>
        parser<S,U> for_each(parser<S,T> p, X k)
        {
            static_assert(std::is_convertible<X, std::function<parser<S,U>(T)>>::value, 
                          "for_each requires a function type parser<U>(T)");
            
            parser<S,U> q = [=](state s) {
                pair_list<S,T> xs = parse(p, s);
                auto yss = functor<list>::fmap<pair_list<S,U>>([=](std::pair<T,S> ps) {
                    T t = ps.first;
                    S s1 = ps.second;
                    parser<S,U> kt = k(t);
                    pair_list<S,U> r = parse(kt, s1);
                    return r;
                }, xs);
                return pair_list<S,U>::flatten(yss);
            };
            return q;
        }

        template <typename T, typename S>
        parser<S,T> fail()
        {
            return [](S s) {
                return pair_list<S,T>{};
            };
        }

        template <typename T, typename S>
        parser<S,T> orelse(parser<S,T> lhs, parser<S,T> rhs)
        {
            return [=](S s) {
                pair_list<S,T> lhs_res = parse(lhs, s);
                if (!lhs_res.empty())
                {
                    return lhs_res;
                }
                return parse(rhs, s);
            };
        }

        template <typename T, typename S, typename X>
        parser<S,T> with(parser<S,T> p, X pred)
        {
            static_assert(std::is_convertible<X, std::function<bool(T)>>::value, 
                          "with requires a function type bool(T)");
            return [=](S s) {
                return filter([=](std::pair<T,S> x) {
                    return pred(x.first);
                }, parse(p, s));
            };
        }

        template <typename S>
        parser<S,char> one()
        {
            return [](S s) {
                if (s.begin() == s.end())
                {
                    return pair_list<S,char>{};
                }
                char r = *s.begin();
                S s1(s.begin() + 1, s.end());
                return pair_list<S,char>(std::make_pair(r, s1));
            };
        }
    }
    // namespace funds
}
// namespace varoom

#endif // VAROOM_FUNDS_PARSING_HPP
