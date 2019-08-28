//
// Created by xiamr on 8/28/19.
//

#include <boost/range/algorithm.hpp>
#include <boost/range/numeric.hpp>
#include <boost/contract.hpp>

#include "CrossCorrelation.hpp"
#include "common.hpp"
#include "Histrogram2D.hpp"

namespace std {
    template<typename T>
    inline std::istream &operator>>(std::istream &is, std::tuple<T, T, T> &t) {
        is >> std::get<0>(t) >> std::get<1>(t) >> std::get<2>(t);
        return is;
    }
}

void CrossCorrelation::calculate(const std::string &out) {
    auto series = readAngleSeries();

    auto correlation_factor = calCovariance(series);

    std::ofstream os(out);
    os << "# " << title() << '\n';

    os << "correlation factor from Covariance Analysis: " << correlation_factor << '\n';

    printCrossCorrelationFunction(calculateCrossCorrelation(series), series, os);

    printConvolutionFunction(calculateConvolutionFunction(series), series, os);

    printHistrogramDistribution(getDistribution(series), os);

    os << std::string(50, '#') << '\n';

}

void
CrossCorrelation::printHistrogramDistribution(const std::vector<std::tuple<double, double, double>> &distribution,
                                              std::ofstream &os) {
    os << std::string(50, '#') << '\n';

    os << boost::format("%15s %15s %15s\n") % "Angle1(degree)" % "Angle2(degree)" % "Probability Density(degree-2)";

    for (auto[angle1, angle2, density] : distribution) {
        os << boost::format("%15.6f %15.5f %15.5f\n") % angle1 % angle2 % (100 * density);
    }

}

void CrossCorrelation::printConvolutionFunction(const std::vector<double> &convolution_function,
                                                std::deque<std::tuple<double, double, double>> &series,
                                                std::ofstream &os) {
    os << std::string(50, '#') << '\n';

    os << boost::format("%15s %20s\n") % "Time" % "convolution function";

    for (std::size_t i = 0; i < convolution_function.size(); ++i) {
        os << boost::format("%15.6f %20.5f\n")
              % (i == 0 ? 0.0 : std::get<0>(series[i - 1]))
              % convolution_function[i];
    }
}

void CrossCorrelation::printCrossCorrelationFunction(const std::vector<double> &cross_correlation_function,
                                                     std::deque<std::tuple<double, double, double>> &series,
                                                     std::ofstream &os) {
    os << std::string(50, '#') << '\n';

    os << boost::format("%15s %20s\n") % "Time" % "Cross-correlation function";

    for (std::size_t i = 0; i < cross_correlation_function.size(); ++i) {
        os << boost::format("%15.6f %20.5f\n")
              % (i == 0 ? 0.0 : std::get<0>(series[i - 1]))
              % cross_correlation_function[i];
    }
}

std::vector<std::tuple<double, double, double>>
CrossCorrelation::getDistribution(const std::deque<std::tuple<double, double, double>> &series) {
    Histrogram2D histrogram;
    histrogram.initialize(180, 5, 180, 5);

    for (auto[time, x, y] : series) {
        histrogram.update(x, y);
    }

    return histrogram.getDistribution();
}

std::deque<std::tuple<double, double, double>> CrossCorrelation::readAngleSeries() {
    auto file = choose_file("Enter Data File: ", true);

    std::ifstream f{file};
    return std::deque<std::tuple<double, double, double>>(
            std::istream_iterator<std::tuple<double, double, double>>{f},
            std::istream_iterator<std::tuple<double, double, double>>{});
}

std::vector<double>
CrossCorrelation::calculateConvolutionFunction(const std::deque<std::tuple<double, double, double>> &series) {
    std::vector<double> convolution_function(series.size());
    std::vector<std::size_t> ntime(series.size());

    for (std::size_t i = 0; i < series.size(); ++i) {
        for (std::size_t t = i; t < series.size(); ++t) {

            convolution_function[t] = std::get<1>(series[i]) * std::get<2>(series[t - i]);
            ++ntime[t];
        }
    }

    for (std::size_t n = 0; n < convolution_function.size(); ++n) {
        convolution_function[n] /= ntime[n];
    }

    return convolution_function;
}

std::vector<double>
CrossCorrelation::calculateCrossCorrelation(std::deque<std::tuple<double, double, double>> &series) {

    std::vector<double> cross_correlation_function(series.size());
    std::vector<std::size_t> ntime(series.size());

    for (std::size_t i = 0; i < series.size(); ++i) {
        for (std::size_t j = i; j < series.size(); ++j) {
            auto n = j - i;
            cross_correlation_function[n] += std::get<1>(series[i]) * std::get<2>(series[j]);
            ++ntime[n];
        }
    }

    for (std::size_t n = 0; n < cross_correlation_function.size(); ++n) {
        cross_correlation_function[n] /= ntime[n];
    }
    return cross_correlation_function;
}

double CrossCorrelation::calCovariance(const std::deque<std::tuple<double, double, double>> &series) {

    boost::contract::check check = boost::contract::public_function<CrossCorrelation>()
            .precondition([&] { BOOST_CONTRACT_ASSERT(!series.empty()); });

    auto[u_sum, v_sum] = boost::accumulate(series, std::pair<double, double>(),
                                           [](auto &p, auto &t) {
                                               std::get<0>(p) += std::get<1>(t);
                                               std::get<1>(p) += std::get<2>(t);
                                               return p;
                                           });
    auto u = u_sum / series.size();
    auto v = v_sum / series.size();

    auto[stdev_u, stdev_v] = boost::accumulate(series, std::pair<double, double>(),
                                               [u, v](auto &p, auto &t) {
                                                   std::get<0>(p) += std::pow(std::get<1>(t) - u, 2);
                                                   std::get<1>(p) += std::pow(std::get<2>(t) - v, 2);
                                                   return p;
                                               });

    stdev_u = std::sqrt(stdev_u / series.size());
    stdev_v = std::sqrt(stdev_v / series.size());


    auto cov = boost::accumulate(series, 0.0, [u, v](auto cov, auto &t) {
        return cov + (std::get<1>(t) - u) * (std::get<2>(t) - v);
    });

    cov /= series.size();

    auto correlation_factor = cov / (stdev_u * stdev_v);

    return correlation_factor;
}
