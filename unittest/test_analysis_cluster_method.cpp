//
// Created by xiamr on 6/20/19.
//

#include <gmock/gmock.h>
#include "../src/Cluster.hpp"


using namespace testing;
using namespace std;


bool operator==(const Cluster::conf_clust &c1, const Cluster::conf_clust &c2) {
    return c1.conf == c2.conf && c1.clust == c2.clust;
}

ostream &operator<<(ostream &out, const Cluster::conf_clust &c) {
    out << "{conf:" << c.conf << " clust:" << c.clust << "}";
    return out;
}


TEST(ClusterTest, initialize_conf_clust_vector) {
    class ClusterTest : public Cluster {
    public:
        std::vector<Cluster::conf_clust> test_initialize_conf_clust_vector(int conf_size) const {
            return Cluster::initialize_conf_clust_vector(conf_size);
        }
    } cluster;

    ASSERT_EQ(cluster.test_initialize_conf_clust_vector(-1).size(), 0);
    ASSERT_EQ(cluster.test_initialize_conf_clust_vector(0).size(), 0);

    ASSERT_THAT(cluster.test_initialize_conf_clust_vector(1), ContainerEq(vector<Cluster::conf_clust>{{0, 0}}));
    ASSERT_THAT(cluster.test_initialize_conf_clust_vector(3),
                ContainerEq(vector<Cluster::conf_clust>{{0, 0},
                                                        {1, 1},
                                                        {2, 2}}));
}

TEST(ClusterTest, do_clust) {
    class ClusterTest : public Cluster {
    public:
        std::vector<Cluster::conf_clust>
        test_do_cluster(const std::list<Cluster::rmsd_matrix> &rmsd_list, int conf_size) const {
            return Cluster::do_cluster(rmsd_list, conf_size);
        }

        explicit ClusterTest(double cutoff) : Cluster() {
            setCutoff(cutoff);
        }
    } cluster(1.0);


    ASSERT_THAT(cluster.test_do_cluster(list<Cluster::rmsd_matrix>{{0, 1, 1.1}}, 2),
                ContainerEq(vector<Cluster::conf_clust>{{0, 0},
                                                        {1, 1}}));

    /*     0     1      2       3
     *  0       1.1    1.2     1.3
     *
     *  1              1.2     1.1
     *
     *  2                      2.0
     *
     */
    ASSERT_THAT(cluster.test_do_cluster(list<Cluster::rmsd_matrix>{{0, 1, 1.1},
                                                                   {1, 3, 1.1},
                                                                   {0, 2, 1.2},
                                                                   {1, 2, 1.2},
                                                                   {0, 3, 1.3},
                                                                   {2, 3, 2.0}}, 4),
                ContainerEq(vector<Cluster::conf_clust>{{0, 0},
                                                        {1, 1},
                                                        {2, 2},
                                                        {3, 3}}));
    /*     0     1      2       3
     *  0       0.8    1.2     1.3
     *
     *  1              1.2     1.1
     *
     *  2                      0.9
     *
     */
    ASSERT_THAT(cluster.test_do_cluster(list<Cluster::rmsd_matrix>{{2, 3, 0.9},
                                                                   {0, 1, 0.8},
                                                                   {1, 3, 1.1},
                                                                   {0, 2, 1.2},
                                                                   {1, 2, 1.2},
                                                                   {0, 3, 1.3}}, 4),
                ContainerEq(vector<Cluster::conf_clust>{{0, 0},
                                                        {1, 0},
                                                        {2, 2},
                                                        {3, 2}}));
}

