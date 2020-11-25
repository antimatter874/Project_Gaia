#include <iterator>
#include <numeric>
#include <vector>

#include "extern/csv.hpp"
#include "impl/maths.hpp"


inline auto prior_instance = [](double x, double par_err, double w_true, int row_count){
    return gaia::maths::gaussian(1/x, par_err, w_true) * (1.0 / row_count);
};

auto pdf(const std::vector<double>& vec,
         const std::function<double(double, double, double, int)>& fn,
         const int n) -> std::vector<std::vector<double>>
{
    static constexpr auto r_true = 100;
    static constexpr auto w_true = 1.0 / r_true;

    auto mat = std::vector<std::vector<double>>(n, std::vector<double>());
    auto col = 0;

    for (auto f : vec)
    {
        auto par_err = f * w_true;

        for (int i = 0; i < n; ++i)
        {
            if (col == 0) { mat.at(i).push_back(i); }
            auto val = fn(i, par_err, w_true, n);
            mat.at(i).push_back(val);
        }
        ++col;
    }

    return mat;
}

void div_max_col(std::vector<std::vector<double>>& mat)
{
    if (mat.empty() or mat[0].empty()) { return; }
    const auto col_count = mat[0].size();

    constexpr auto min = std::numeric_limits<double>::min();
    auto max_col_values = std::vector<double>(col_count - 1, min);

    for (auto const& row : mat)
    {
        for (std::size_t i = 1; i < col_count; ++i)
        {
            auto col = row.at(i);
            auto& max_val = max_col_values[i];
            if (col > max_val) { max_val = col; }
        }
    }

    for (auto& row : mat)
    {
        for (std::size_t i = 1; i < col_count; ++i)
        {
            auto& col = row[i];
            auto max_val = max_col_values[i];
            col /= max_val;
        }
    }
}

int main()
{
    // setup vector
    auto frc_errs = std::vector<double>();
    auto err_count = 6;
    frc_errs.reserve(err_count);
    for (int i = 1; i <= err_count; ++i) { frc_errs.push_back(i / 10.0); }

    // generate matrix
    auto row_count = 600;
    auto mat = pdf(frc_errs, prior_instance, row_count);
    div_max_col(mat);

    // setup csv
    auto path = "dg.csv";
    auto ofile = std::ofstream(path);
    auto writer = csv::make_csv_writer(ofile);

    // write mat to csv
    for (auto& vec : mat) { writer << vec; }
}
