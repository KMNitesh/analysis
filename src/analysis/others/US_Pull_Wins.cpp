
#include "others/US_Pull_Wins.hpp"
#include "utils/common.hpp"

void US_Pull_Wins::process() {
    std::string filename = choose_file("Choose <summary_distances.dat> file : ")
                               .can_empty(true)
                               .extension("dat")
                               .isExist(true)
                               .default_("summary_distances.dat");

    struct Item {
        Item(int conf_number, double distance) : conf_number(conf_number), distance(distance) {}

        int conf_number;
        double distance;
    };
    std::vector<Item> item_vector;

    std::ifstream ifs(filename);
    std::string line;
    while (!ifs.eof()) {
        std::getline(ifs, line);

        auto fields = split(line);
        if (fields.empty())
            break;
        item_vector.emplace_back(std::stoi(fields[0]), std::stod(fields[1]));
    }

    std::stable_sort(begin(item_vector), end(item_vector),
                     [](const Item &lhs, const Item &rhs) { return lhs.distance < rhs.distance; });

    std::cout << "Begin : " << item_vector.front().distance << '\n';
    std::cout << "  End : " << item_vector.back().distance << '\n';

    double step_size = choose(0.01, "Step size [0.05 nm] : ", Default(0.05));
    int total_windows = choose(1, "Windows number : ");

    double obj_distance = item_vector.front().distance;

    std::vector<Item> find_;

    for (auto it = begin(item_vector); total_windows-- and it != std::end(item_vector); obj_distance += step_size) {
        for (; it != std::end(item_vector);) {
            auto delta = it->distance - obj_distance;
            auto delta_abs = std::abs(delta);
            if (delta_abs < 0.001) {
                find_.push_back(*it);
                break;
            }
            if (delta < 0) {
                ++it;
            } else {
                auto prev_delta_abs = std::abs((it - 1)->distance - obj_distance);
                find_.push_back(delta_abs < prev_delta_abs ? *it : *(it - 1));
                break;
            }
        }
    }

    std::cout << "Total Finding :\n";
    for (const auto &i : find_) {
        std::cout << i.conf_number << "  " << i.distance << '\n';
    }
    std::cout << "List > ";
        for (const auto &i : find_) {
        std::cout << " " << i.conf_number;
    }
    std::cout << std::endl;
}
