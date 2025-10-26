#include <algorithm>
#include <bit>
#include <cmath>
#include <filesystem>
#include <format>
#include <fstream>
#include <future>
#include <limits>
#include <print>
#include <queue>
#include <random>
#include <set>
#include <stdint.h>
#include <thread>

#include "cubic_bond_percolation.h"

#include "colour_names.h"
#include "gnuplot-iostream.h"
#include "pcg_extras.hpp"
#include "pcg_random.hpp"
#include "percolation.h"
#include "power.h"
#include "timer.h"

// GNU plot has its limitations here. Do not waste too much time fiddling with it, will probably write something proper later anyway.

cubic_bond_percolation::cubic_bond_percolation(uint8_t cube_pow, double p)
    : percolation(ipow(2, cube_pow * 3u)), _cube_pow(cube_pow), _cube_size(ipow(2, cube_pow)), _probability(p),
      _bound(std::numeric_limits<uint64_t>::max() * p), _rng(pcg_extras::seed_seq_from<std::random_device>{})
{
  _gp << "set xrange [0:" << _cube_size << "]" << std::endl;
  _gp << "set yrange [0:" << _cube_size << "]" << std::endl;
  _gp << "set zrange [0:" << _cube_size << "]" << std::endl;
  _gp << "set view equal xyz" << std::endl; // Force grid to be square
  _gp << "unset border" << std::endl;
  _gp << "unset xtics" << std::endl;
  _gp << "unset ytics" << std::endl;
  _gp << "unset ztics" << std::endl;
  _gp << "set key outside right top samplen 2 spacing .7 font ',8' tc rgb 'grey40'" << std::endl;
}

void cubic_bond_percolation::set_probability(double p)
{
  _probability = p;
  _bound = std::numeric_limits<uint64_t>::max() * p;
}

// For now, we only use 2^n threads, and max_num_threads is assumed to be >= 2
void cubic_bond_percolation::generate_clusters_parallel(uint8_t max_num_threads)
{
  generate_merge_clusters_recursive(max_num_threads, 0, _cube_size);

  return;
}

void cubic_bond_percolation::generate_merge_clusters_recursive(uint8_t max_num_threads, int start_i, int end_i)
{
  int middle_i = (start_i + end_i) / 2;

  if (std::min(middle_i - start_i, end_i - middle_i) >= 2 * _cube_size / max_num_threads)
  {
    // Split each in half again and recurse, joining up afterwards
    std::thread t1(&cubic_bond_percolation::generate_merge_clusters_recursive, this, max_num_threads, start_i, middle_i);
    std::thread t2(&cubic_bond_percolation::generate_merge_clusters_recursive, this, max_num_threads, middle_i, end_i);

    t1.join();
    t2.join();

    merge_clusters_slices(middle_i);
    return;
  }

  std::thread t1(&cubic_bond_percolation::generate_clusters_parallel_thread, this, start_i, middle_i);
  std::thread t2(&cubic_bond_percolation::generate_clusters_parallel_thread, this, middle_i, end_i);

  t1.join();
  t2.join();

  merge_clusters_slices(middle_i);

  return;
}

void cubic_bond_percolation::generate_clusters_parallel_thread(int start_i, int end_i)
{
  // std::println("Generating clusters for i {}-{}", start_i, end_i);

  pcg64_fast rng(pcg_extras::seed_seq_from<std::random_device>{});

  std::tuple<int, int, int> new_node;

  for (std::get<0>(new_node) = start_i; std::get<0>(new_node) < end_i; ++std::get<0>(new_node))
  {
    for (std::get<1>(new_node) = 0; std::get<1>(new_node) < _cube_size; ++std::get<1>(new_node))
    {
      for (std::get<2>(new_node) = 0; std::get<2>(new_node) < _cube_size; ++std::get<2>(new_node))
      {
        // Make set and merge with previous (up to three) if there is an edge between them
        this->make_set(new_node);

        // Do not make this a loop as no need to construct nodes if rng() >= _bound
        if (rng() < _bound)
        {
          cubic_bond_percolation::merge({std::get<0>(new_node), std::get<1>(new_node), std::max(std::get<2>(new_node) - 1, 0)}, new_node);
        }
        if (rng() < _bound)
        {
          cubic_bond_percolation::merge({std::get<0>(new_node), std::max(std::get<1>(new_node) - 1, 0), std::get<2>(new_node)}, new_node);
        }
        if (rng() < _bound)
        {
          cubic_bond_percolation::merge({std::max(std::get<0>(new_node) - 1, start_i), std::get<1>(new_node), std::get<2>(new_node)}, new_node);
        }
      }
    }
  }

  // std::println("Finished generating clusters for i {}-{}", start_i, end_i);

  return;
}

void cubic_bond_percolation::merge_clusters_slices(int i)
{
  pcg64_fast rng(pcg_extras::seed_seq_from<std::random_device>{});

  // std::println("Merging clusters for i={}", i);
  for (int j = 0; j < _cube_size; ++j)
  {
    for (int k = 0; k < _cube_size; ++k)
    {
      // Make set and merge with previous (up to three) if there is an edge between them
      std::tuple<int, int, int> node1(i, j, k);
      std::tuple<int, int, int> node2(i - 1, j, k);

      if (rng() < _bound)
      {
        cubic_bond_percolation::merge(node1, node2);
      }
    }
  }

  // std::println("Finished merging clusters for i={}", i);

  return;
}

void cubic_bond_percolation::generate_clusters()
{
  // std::println("Generating clusters...");

  std::tuple<int, int, int> new_node;

  std::array<std::tuple<int, int, int>, 3> previous_nodes;

  for (std::get<0>(new_node) = 0; std::get<0>(new_node) < _cube_size; ++std::get<0>(new_node))
  {
    for (std::get<1>(new_node) = 0; std::get<1>(new_node) < _cube_size; ++std::get<1>(new_node))
    {
      for (std::get<2>(new_node) = 0; std::get<2>(new_node) < _cube_size; ++std::get<2>(new_node))
      {
        // Make set and merge with previous (up to three) if there is an edge between them
        this->make_set(new_node);

        previous_nodes[0] = {std::get<0>(new_node), std::get<1>(new_node), std::max(std::get<2>(new_node) - 1, 0)};
        previous_nodes[1] = {std::get<0>(new_node), std::max(std::get<1>(new_node) - 1, 0), std::get<2>(new_node)};
        previous_nodes[2] = {std::max(std::get<0>(new_node) - 1, 0), std::get<1>(new_node), std::get<2>(new_node)};

        for (const auto& node : previous_nodes)
        {
          if (_rng() < _bound)
          {
            cubic_bond_percolation::merge(node, new_node);
          }
        }
      }
    }
  }

  // std::println("Finished generating clusters");

  return;
}

void cubic_bond_percolation::plot_clusters(uint32_t min_cluster_size, size_t max_num_clusters) const
{
  max_num_clusters = std::min(colour_names.size(), max_num_clusters);
  _gp << "set title tc rgb 'grey40' 'Percolation, p=" << std::setprecision(8) << _probability << " Cube size=" << _cube_size << "'" << std::endl;

  const auto largest_clusters = this->get_clusters_sorted(min_cluster_size);
  size_t count = 0;

  // std::println("Number of clusters of size at least {}: {}", min_cluster_size, largest_clusters.size());

  for (auto it = largest_clusters.crbegin(); it != largest_clusters.crend() && count < max_num_clusters; ++it, ++count)
  {
    _gp << ((count == 0) ? "splot" : "replot") << _gp.file1d(it->second) << "u 1:2:3:(0.03) with points lc rgb '" << colour_names[count]
        << "' pt 7 ps variable title 'Cluster " << count + 1 << " (" << it->second.size() << " points)"
        << ((it->first.size > 0) ? "(terminated)" : "(still growing)") << "'"
        << ((count == max_num_clusters - 1 || it == --largest_clusters.crend()) ? "; pause mouse close" : "") << std::endl;
  }
}

void cubic_bond_percolation::plot_central_clusters(uint32_t min_cluster_size, size_t central_cube_size, size_t max_num_clusters) const
{
  max_num_clusters = std::min(colour_names.size(), max_num_clusters);

  std::tuple<int, int, int> current_node;
  std::map<node, std::vector<std::tuple<int, int, int>>> clusters;

  // Populate map with all roots of clusters which intersect a central cube
  if (central_cube_size > _cube_size)
  {
    std::println("Central cube larger than simulation");
    return;
  }

  const size_t min_coordinate = (_cube_size - central_cube_size) / 2;
  const size_t max_coordinate = (_cube_size + central_cube_size) / 2;
  for (std::get<0>(current_node) = min_coordinate; std::get<0>(current_node) < max_coordinate; ++std::get<0>(current_node))
  {
    for (std::get<1>(current_node) = min_coordinate; std::get<1>(current_node) < max_coordinate; ++std::get<1>(current_node))
    {
      for (std::get<2>(current_node) = min_coordinate; std::get<2>(current_node) < max_coordinate; ++std::get<2>(current_node))
      {
        const size_t index = this->get_index(current_node);

        const node& root = *this->find_const(&this->_forest[index]);

        if (std::abs(root.size) >= min_cluster_size)
        {
          if (!clusters.contains(root))
          {
            clusters.emplace(root, std::vector<std::tuple<int, int, int>>({this->get_element(index)}));
          }
        }
      }
    }
  }

  // Now populate map fully
  for (std::get<0>(current_node) = 0; std::get<0>(current_node) < _cube_size; ++std::get<0>(current_node))
  {
    for (std::get<1>(current_node) = 0; std::get<1>(current_node) < _cube_size; ++std::get<1>(current_node))
    {
      for (std::get<2>(current_node) = 0; std::get<2>(current_node) < _cube_size; ++std::get<2>(current_node))
      {
        const size_t index = this->get_index(current_node);

        const node& root = *this->find_const(&this->_forest[index]);

        if (std::abs(root.size) >= min_cluster_size)
        {
          if (clusters.contains(root))
          {
            clusters.at(root).push_back(this->get_element(index));
          }
        }
      }
    }
  }

  _gp << "set title tc rgb 'grey40' 'Percolation, p=" << std::setprecision(8) << _probability << " Cube size=" << _cube_size << "'" << std::endl;

  size_t count = 0;

  // std::println("Number of clusters of size at least {}: {}", min_cluster_size, clusters.size());

  for (auto it = clusters.crbegin(); it != clusters.crend() && count < max_num_clusters; ++it, ++count)
  {
    _gp << ((count == 0) ? "splot" : "replot") << _gp.file1d(it->second) << "u 1:2:3:(0.03) with points lc rgb '" << colour_names[count]
        << "' pt 7 ps variable title 'Cluster " << count + 1 << " (" << it->second.size() << " points)"
        << ((it->first.size > 0) ? "(terminated)" : "(still growing)") << "'"
        << ((count == max_num_clusters - 1 || it == --clusters.crend()) ? "; pause mouse close" : "") << std::endl;
  }
}

void cubic_bond_percolation::write_clusters_data(uint32_t min_cluster_size, size_t central_cube_size) const
{
  std::tuple<int, int, int> current_node;
  std::set<node> clusters;

  // Populate map with all roots of clusters which intersect a central cube
  if (central_cube_size > _cube_size)
  {
    std::println("Central cube larger than simulation");
    return;
  }

  const size_t min_coordinate = (_cube_size - central_cube_size) / 2;
  const size_t max_coordinate = (_cube_size + central_cube_size) / 2;
  for (std::get<0>(current_node) = min_coordinate; std::get<0>(current_node) < max_coordinate; ++std::get<0>(current_node))
  {
    for (std::get<1>(current_node) = min_coordinate; std::get<1>(current_node) < max_coordinate; ++std::get<1>(current_node))
    {
      for (std::get<2>(current_node) = min_coordinate; std::get<2>(current_node) < max_coordinate; ++std::get<2>(current_node))
      {
        const size_t index = this->get_index(current_node);

        const node& root = *this->find_const(&this->_forest[index]);

        if (std::abs(root.size) >= min_cluster_size)
        {
          if (!clusters.contains(root))
          {
            clusters.emplace(root);
          }
        }
      }
    }
  }

  // TODO: ensure directory exists
  std::ofstream data_file(
      std::format("src/analyse_data/data/test/cubic_bond_percolation_p_{:.10f}_centre_{}_size_{}.csv", _probability, central_cube_size, _cube_size));

  data_file << "probability, central cube size, simulation size, number of simulations\n";
  data_file << std::format("{:.10f}, {}, {}, 1\n", _probability, central_cube_size, _cube_size);
  data_file << "\nsize,number terminated,number still growing\n";

  std::array<size_t, 3> line = {static_cast<size_t>(std::abs(clusters.crbegin()->size)), 0, 0};

  for (auto it = clusters.crbegin(); it != clusters.crend(); ++it)
  {
    if (std::abs(it->size) == line[0])
    {
      if (it->size > 0)
      {
        ++line[1];
      }
      else
      {
        ++line[2];
      }
      continue;
    }

    // Finished with last line, so write out
    data_file << std::format("{},{},{}\n", line[0], line[1], line[2]);

    line[0] = std::abs(it->size);
    line[1] = 0;
    line[2] = 0;
    if (it->size > 0)
    {
      ++line[1];
    }
    else
    {
      ++line[2];
    }
  }

  data_file << std::format("{},{},{}\n", line[0], line[1], line[2]);
}

void cubic_bond_percolation::run_simulations(uint32_t num_simulations, size_t central_cube_size)
{
  std::println("Running {} simulations with size {} for p={}", num_simulations, _cube_size, _probability);
  timer tm;

  if (central_cube_size > _cube_size)
  {
    std::println("Central cube larger than simulation");
    return;
  }

  std::vector<std::pair<uint64_t, uint64_t>> results;
  for (uint32_t simulation_count = 0; simulation_count < num_simulations; ++simulation_count)
  {
    std::print("Simulation number: {}", simulation_count);
    tm.restart();
    // Run simulation
    generate_clusters_parallel(4);

    std::tuple<int, int, int> current_node;
    std::set<node> clusters; // Should really make this unordered

    // Populate map with all roots of clusters which intersect a central cube
    const size_t min_coordinate = (_cube_size - central_cube_size) / 2;
    const size_t max_coordinate = (_cube_size + central_cube_size) / 2;
    for (std::get<0>(current_node) = min_coordinate; std::get<0>(current_node) < max_coordinate; ++std::get<0>(current_node))
    {
      for (std::get<1>(current_node) = min_coordinate; std::get<1>(current_node) < max_coordinate; ++std::get<1>(current_node))
      {
        for (std::get<2>(current_node) = min_coordinate; std::get<2>(current_node) < max_coordinate; ++std::get<2>(current_node))
        {
          const size_t index = this->get_index(current_node);

          const node& root = *this->find_const(&this->_forest[index]);

          if (!clusters.contains(root))
          {
            clusters.emplace(root);
          }
        }
      }
    }

    // Ensure we have space to store all results
    uint32_t largest_bucket = std::bit_width(static_cast<uint32_t>(std::abs(clusters.crbegin()->size))) - 1;
    if (results.size() < largest_bucket + 1)
    {
      results.resize(largest_bucket + 1, std::pair<uint64_t, uint64_t>(0, 0));
    }

    // Populate results
    uint64_t total_size = 0;
    for (auto it = clusters.cbegin(); it != clusters.cend(); ++it)
    {
      /*
      We might like to be able to have more data points by using a smaller base for log, however this leads to discrepancies
      as the different buckets no longer have integer boundaries. So we use log base 2.
      */
      uint32_t bucket = std::bit_width(static_cast<uint32_t>(std::abs(it->size))) - 1;
      if (it->size > 0)
      {
        ++results[bucket].first;
        total_size += it->size;
      }
      else
      {
        ++results[bucket].second;
        total_size += -it->size;
      }
    }

    tm.stop();
    std::print(" finished in {} ms\n", tm.get_ms());
  }

  // Write out results to file
  std::filesystem::path results_path = std::format("src/analyse_data/data/p_244_test2/cubic_bond_percolation_p_{:.10f}_centre_{}_size_{}_num_{}.csv",
                                                   _probability, central_cube_size, _cube_size, num_simulations);
  std::filesystem::create_directory(results_path.parent_path());
  std::ofstream data_file(results_path);

  data_file << "probability, central cube size, simulation size, number of simulations\n";
  data_file << std::format("{:.10f}, {}, {}, {}\n", _probability, central_cube_size, _cube_size, num_simulations);
  data_file << "\nstart size,number terminated,number still growing\n";

  for (size_t bucket = 0; bucket < results.size(); ++bucket)
  {
    data_file << std::format("{}, {}, {}\n", bucket + 1, results[bucket].first, results[bucket].second);
  }

  std::println("Completed {} simulations with size {} for p={}", num_simulations, _cube_size, _probability);
}

void cubic_bond_percolation::run_simulations_test(const std::string& folder_name, uint32_t num_simulations, size_t central_cube_size)
{
  std::println("Running {} simulations with size {} for p={}", num_simulations, _cube_size, _probability);
  timer tm;

  if (central_cube_size > _cube_size)
  {
    std::println("Central cube larger than simulation");
    return;
  }

  std::vector<std::pair<uint64_t, uint64_t>> results;
  for (uint32_t simulation_count = 0; simulation_count < num_simulations; ++simulation_count)
  {
    std::print("Simulation number: {}", simulation_count);
    tm.restart();
    // Run simulation
    generate_clusters_parallel(4);

    auto new_results =
        count_clusters_parallel_recursive(4, (_cube_size - central_cube_size) / 2, (_cube_size + central_cube_size) / 2, central_cube_size);

    if (new_results.size() > results.size())
    {
      results.resize(new_results.size(), std::pair<uint64_t, uint64_t>(0, 0));
    }

    for (size_t i = 0; i < new_results.size(); ++i)
    {
      results[i].first += new_results[i].first;
      results[i].second += new_results[i].second;
    }

    tm.stop();
    std::print(" finished in {} ms\n", tm.get_ms());
  }

  // Write out results to file
  std::filesystem::path results_path = std::format("src/analyse_data/data/{}/cubic_bond_percolation_p_{:.10f}_centre_{}_size_{}_num_{}.csv",
                                                   folder_name, _probability, central_cube_size, _cube_size, num_simulations);
  std::filesystem::create_directory(results_path.parent_path());
  std::ofstream data_file(results_path);

  data_file << "probability, central cube size, simulation size, number of simulations\n";
  data_file << std::format("{:.10f}, {}, {}, {}\n", _probability, central_cube_size, _cube_size, num_simulations);
  data_file << "\nstart size,number terminated,number still growing\n";

  for (size_t bucket = 0; bucket < results.size(); ++bucket)
  {
    data_file << std::format("{}, {}, {}\n", bucket + 1, results[bucket].first, results[bucket].second);
  }

  std::println("Completed {} simulations with size {} for p={}", num_simulations, _cube_size, _probability);
}

std::vector<std::pair<uint64_t, uint64_t>> cubic_bond_percolation::count_clusters_parallel_recursive(uint8_t max_num_threads, int start_i, int end_i,
                                                                                                     size_t central_cube_size) const
{
  std::vector<std::pair<uint64_t, uint64_t>> results1;
  std::vector<std::pair<uint64_t, uint64_t>> results2;

  int middle_i = (start_i + end_i) / 2;

  if (std::min(middle_i - start_i, end_i - middle_i) >= 2 * central_cube_size / max_num_threads)
  {

    // Split each in half again and recurse, joining up afterwards
    std::future<std::vector<std::pair<uint64_t, uint64_t>>> promise1 = std::async(
        std::launch::async, &cubic_bond_percolation::count_clusters_parallel_recursive, this, max_num_threads, start_i, middle_i, central_cube_size);
    std::future<std::vector<std::pair<uint64_t, uint64_t>>> promise2 = std::async(
        std::launch::async, &cubic_bond_percolation::count_clusters_parallel_recursive, this, max_num_threads, middle_i, end_i, central_cube_size);

    results1 = promise1.get();
    results2 = promise2.get();
  }
  else
  {
    std::future<std::vector<std::pair<uint64_t, uint64_t>>> promise1 =
        std::async(std::launch::async, &cubic_bond_percolation::count_clusters_parallel_thread, this, start_i, middle_i, central_cube_size);
    std::future<std::vector<std::pair<uint64_t, uint64_t>>> promise2 =
        std::async(std::launch::async, &cubic_bond_percolation::count_clusters_parallel_thread, this, middle_i, end_i, central_cube_size);

    results1 = promise1.get();
    results2 = promise2.get();
  }

  // Merge results and return
  if (results2.size() > results1.size())
  {
    results1.resize(results2.size(), std::pair<uint64_t, uint64_t>(0, 0));
  }

  for (size_t i = 0; i < results2.size(); ++i)
  {
    results1[i].first += results2[i].first;
    results1[i].second += results2[i].second;
  }

  return results1;
}

std::vector<std::pair<uint64_t, uint64_t>> cubic_bond_percolation::count_clusters_parallel_thread(int start_i, int end_i,
                                                                                                  size_t central_cube_size) const
{
  std::vector<std::pair<uint64_t, uint64_t>> results;

  std::tuple<int, int, int> current_node;
  // Populate map with all roots of clusters which intersect a central cube
  const size_t min_coordinate = (_cube_size - central_cube_size) / 2;
  const size_t max_coordinate = (_cube_size + central_cube_size) / 2;
  for (std::get<0>(current_node) = start_i; std::get<0>(current_node) < end_i; ++std::get<0>(current_node))
  {
    for (std::get<1>(current_node) = min_coordinate; std::get<1>(current_node) < max_coordinate; ++std::get<1>(current_node))
    {
      for (std::get<2>(current_node) = min_coordinate; std::get<2>(current_node) < max_coordinate; ++std::get<2>(current_node))
      {
        const size_t index = this->get_index(current_node);

        const node& root = *this->find_const(&this->_forest[index]);

        uint32_t bucket = std::bit_width(static_cast<uint32_t>(std::abs(root.size))) - 1;
        if (results.size() < bucket + 1)
        {
          results.resize(bucket + 1, std::pair<uint64_t, uint64_t>(0, 0));
        }

        if (root.size > 0)
        {
          ++results[bucket].first;
        }
        else
        {
          ++results[bucket].second;
        }
      }
    }
  }

  return results;
}

int main()
{
  cubic_bond_percolation perc(9, 0.24); // p(0.2488125); 2^8 = 256 2^9 = 512

  // TODO: plots show size as being one too large.

  for (auto [probability, count] = std::tuple<double, size_t>{0.24878, 0}; count < 1; ++count, probability += 0.00001)
  {
    std::println("Loop {}: Generating clusters for probability={:.10f}", count, probability);
    perc.set_probability(probability);
    perc.run_simulations_test("p_244_test25", 20, 32);

    // perc.generate_clusters_parallel(4);
    // perc.write_clusters_data(1, 64);
  }
  // perc.plot_central_clusters(1000, 16, 10);
  return 0;
}

/*
If we want to use mmap, it is much too slow to use directly due to the random access nature of the disjoint set forest.
Perhaps it would be possible to run the simulation for each slice in memory. Then we could mmap these vectors to chunks of a
larger file, before finally merging these chunks to get the final result - massively reducing thrashing
*/
