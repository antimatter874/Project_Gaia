#ifndef DG_ALGORITHMS_HPP
#define DG_ALGORITHMS_HPP


namespace gaia::algos {

    template <class InputIt, class OutputIt>
    OutputIt distance(InputIt begin_par, const InputIt end_par, OutputIt out)
    {
        while (begin_par != end_par) { *out++ = 1 / *begin_par++; }
        return out;
    }

    template <class InputIt, class OutputIt>
    OutputIt distance_error(InputIt begin_par, const InputIt end_par,
                            InputIt begin_par_err, OutputIt out)
    {
        while (begin_par != end_par)
        {
            const auto par_val = *begin_par++;
            *out++ = *begin_par_err++ / (par_val * par_val);
        }
        return out;
    }

    template <class InputIt, class OutputIt>
    OutputIt fractional_error(InputIt begin_par, const InputIt end_par,
                              InputIt begin_par_err, OutputIt out)
    {
        while (begin_par != end_par) { *out++ = *begin_par_err / *begin_par++; }
        return out;
    }

}

#endif //DG_ALGORITHMS_HPP
