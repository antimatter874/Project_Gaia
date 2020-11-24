#include <iostream>
#include <fstream>
#include <string_view>
#include <string>
#include <algorithm>
#include <cstring>
#include <cmath>
#include <sstream>
#include <cstdio>
#include <iostream>
#include <vector>
#include <numeric>
#include <iterator>
#include <functional>
#include <cwchar>

#include "csv.hpp"

using sviews_t = std::vector<std::string_view>;

bool column_checker(std::string_view input_file, const sviews_t& columns_wanted)
{
    auto reader = csv::CSVReader(input_file);
    auto col_names = reader.get_col_names();
    for (const auto& col : columns_wanted)
    {
        if (std::find(col_names.cbegin(), col_names.cend(), col) == col_names.end())
        {
            return false;
        }
    }
}

void output_data(std::string_view input_file, std::string_view output_file, const sviews_t& columns_wanted)
{
    // check columns exist
//    column_checker(input_file, columns_wanted);

    // setup
    auto reader = csv::CSVReader(input_file);
    auto ofile  = std::ofstream(output_file.data());

    // parse csv
    for (auto& row : reader)
    {
        for (const auto& col : columns_wanted) { ofile << row[col.data()].get_sv() << ' '; }
        ofile << '\n';
    }
}

//std:vector load_data(std::string_view input_file, const std::string_view column_wanted)
//{
//    // check columns exist
//    column_checker(input_file, columns_wanted);
//
//    // setup
//    auto reader = csv::CSVReader(input_file);
//    std::vector<float> v;
//
//    // parse csv
//    for (auto& row: reader) {
//        // Note: Can also use index of column with [] operator
//         v.push_back(row[column_wanted].get_sv());
//    }
//
//    return v
//}


std::vector<double> calc_distance(const std::vector<double>& par)
{
    auto ret = std::vector<double>();
    std::transform(par.cbegin(), par.cend(),
                   std::back_inserter(ret),
                   [](double d){ return 1.0 / d;});
    return ret;
}


std::vector<double> calc_distance_error(const std::vector<double>& par,
                                        const std::vector<double>& par_err)
{
    assert(par.size() == par_err.size());
    auto ret = std::vector<double>();
    ret.reserve(par.size());
    for (std::size_t i = 0; i < par.size(); ++i) {
        double dist_err = par_err[i] / (par[i] * par[i]);
        ret.push_back(dist_err);
    }
    return ret;
}

std::vector<double> calc_fractional_err(const std::vector<double>& par,
                                        const std::vector<double>& par_err)
{
    assert(par.size() == par_err.size());
    auto ret = std::vector<double>();
    ret.reserve(par.size());
    for (std::size_t i = 0; i < par.size(); ++i) {
        double frac_err = par_err[i] / (par[i]);
        ret.push_back(frac_err);
    }
    return ret;
}

double gaussian(const double x, const double sigma, const double mu)
{
    return (1.0/(sqrt(2.0*M_PI) * sigma)) * exp(((-1.0) / (2.0 * sigma * sigma)) * ((x - mu) * (x - mu)));
}



double trapzd(const std::function<double(double)>& fn, double a, double b, int n)
{
    if (n == 1)
    {
        double s = 0.5 * (b-a) * fn(a) + fn(b);
        return s;
    }
    else
    {
        int i = (n - 1);  // same as the loop

        double h = (b - a ) / n;
        double x = a;
        double sum = 0.0;

        for (int j = 0; j <= i; ++j) { sum += 2*fn(x); x += h; }

        double s = 0.5 * h * sum;
        return s;
    }
}

double prior(double x)
{
    return double (1.0/(2.0*1000.0*1000.0*1000.0))*x*x*exp(-x/1000.0);
}


void export_data_simple(std::string filename, const std::vector<double>& data)
{
    filename = std::string ("./" + filename);
    std::ofstream output_file(filename);
    for (const auto &i : data) output_file << i << "\n";
}

void export_data(std::string filename, const std::vector<std::vector<double>>& data)
{
    filename = std::string ("./" + filename);
//    auto ofile  = std::ofstream(filename.data());
//    std::ofstream output_file(filename); // Can also use ofstream, etc.
//    auto writer = csv::make_csv_writer(ofile);
//
//    for (auto& v : data) { writer << v; }

    std::ofstream out(filename);

    for (auto& row : data) {
        for (auto col : row)
            out << col <<',';
        out << '\n';
    }

//
//    for (auto& v : data)
//    {
//        for (auto& i : v){output_file << i << " ";}
//        output_file << "\n";
//    };
}

std::vector<double> largestInColumn(const std::vector<std::vector<double>>& data, int rows, int cols)
{
    std::vector<double> max;
    for (int i = 0; i < cols; i++) {
        // initialize the maximum element
        // with 0
        double maxm = data[0][i];

        // Run the inner loop for rows
        for (int j = 1; j < rows; j++) {
            // check if any element is greater
            // than the maximum element
            // of the column and replace it
            if (data[j][i] > maxm)
                maxm = data[j][i];
        }

        // print the largest element of the column
        max.push_back(maxm);
    }

    return max;
}


int main(){


    auto uniform = [](double x) {
        return double(1.0/600);
    };

//    auto uniform_sd = [](double x) {
//        return (3.0/(600*600*600))*(x*x);
//    };

// calculating pdf values for different fractional errors
    std::vector<double> frc_errs{ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6 };

    int n = 600;
    int m = frc_errs.size();

    std::vector<std::vector<double>> prior1_norm;
//    prior1_norm.reserve(n);

    int col = 0;
    for (double f : frc_errs){

        double frac_parr_error1 = f;
        double r_true1 = 100;
        double w_true1 = 1.0 / r_true1;
        double parerr1 = frac_parr_error1 * w_true1;

        auto prior_instance = [parerr1, w_true1, uniform](double x) {
            return gaussian(1/x, parerr1, w_true1)*uniform(x);
        };

//        double area1 = trapzd(prior_instance, 0, 10000, 100000);
        std::cout << col;
        for (int i = 1; i < n; i++) {
            auto r = double(i);

//            if (col == 0){prior1_norm[i].push_back(r);}

            double prior_value = prior_instance(r);
            prior1_norm[i].push_back(prior_value);
        }
        col++;




//        std::vector<double> max = largestInColumn(prior1_norm, 10000, 2);
//        std::cout << max[1];
//        for (int j = 1; j<10000; j++){
//            prior1_norm[j][1] = prior1_norm[j][1]/max;
//        }


    }

//    std::cout << prior1_norm[300][0];
//std::string filename1 = "uniform_norm" + std::to_string(int(f*10)) + ".csv";
    std::string filename1 = "uniform_norm.csv";
    export_data(filename1, prior1_norm);

    prior1_norm.clear();




//    std::vector<std::vector<double>> prior2_norm(10000, std::vector<double>(2));
////        prior1_norm.reserve(10000);
//    double frac_parr_error1 = 0.1;
////        std::cout << f;
//    double r_true1 = 100;
//    double w_true1 = 1.0 / r_true1;
//    double parerr1 = frac_parr_error1 * w_true1;
//
//    auto prior_instance = [parerr1, w_true1](double x) {
//        return gaussian(1/x, parerr1, w_true1)*prior(x);
//    };
//
//    for (int i = 1; i < 10000; i++) {
//        double r = double(i);
//        std::cout << r << "\n";
//        double prior_value = prior_instance(r);
//        prior2_norm.push_back({r, prior_value});
//    }
//
//    export_data("./testing.csv", prior2_norm);

//    std::cout << prior1_norm[100][0];
//    for (int i = 0; i < 10000; i++) {
//        for (int j = 0; j < 2; j++){
//            std::cout << double (prior2_norm[i][j])<< " ";
//            }
//        std::cout<< "\n";
//        }
//    std::cout << prior2_norm[2][100];
}
