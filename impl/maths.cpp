#ifndef DG_MATHS_HPP
#define DG_MATHS_HPP

#include <cassert>
#include <cmath>
#include <functional>
// #include <numbers>

namespace gaia::maths {

    template <typename Real = double>
    inline Real gaussian(Real x, Real sigma, Real mu) noexcept
    {
        // constexpr auto tau = std::numbers::pi_v<Real> * 2;
        constexpr auto tau = static_cast<Real>(M_PI) * 2;
        Real a = 1 / (std::sqrt(tau) * sigma);
        Real b = std::exp((-1 / (2 * sigma * sigma)) * ((x - mu) * (x - mu)));
        return a * b;
    }

    template <typename Real = double, typename Integer = int>
    inline Real trapzd(const std::function<Real(Real)>& fn, Real a, Real b, Integer n) noexcept
    {
        assert(n > 0);
        if (n == 1) { return 0.5 * (b - a) * fn(a) * fn(b); }
        else
        {
            Integer i = n - 1;

            Real h = (b - a) / n;
            Real x = a;
            Real sum = 0;

            for (Integer j = 0; j <= i; ++j) { sum += 2 * fn(x); x += h; }
            return 0.5 * h * sum;
        }
    }

    template <typename Real = double>
    inline Real prior(Real x) noexcept
    {
        Real ret = x * x * std::exp(-x / 1'000);
        return ret / ((Real)2 * 1'000 * 1'000 * 1'000);
    }
}

#endif //DG_MATHS_HPP
