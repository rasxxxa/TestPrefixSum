#include <iostream>
#include <array>
#include <random>
#include <fstream>
#include <map>
#include <functional>
#include <filesystem>
#include <vector>
#include <set>
#include <fstream>
#include <chrono>
#include <unordered_set>
#include <deque>

constexpr size_t M = 7;
constexpr size_t N = 5;


constexpr size_t NUMBER_OF_RECTANGLES = M * N;

struct Rectangle
{
	size_t m_M;
	size_t m_N;
	size_t m_area;
	size_t m_value;
	constexpr Rectangle(size_t m , size_t n, size_t value) : m_M(m), m_N(n), m_area(m_M * m_N), m_value(value)
	{

	}

	friend std::ostream& operator<< (std::ostream& of, const Rectangle& rectangle)
	{
		of << rectangle.m_N << " " << rectangle.m_M << " " << rectangle.m_area << " " << rectangle.m_value << std::endl;
		return of;
	}
};

struct Cluster
{
	size_t i;
	size_t j;
	size_t k;
	size_t l;
	size_t value;
	Cluster(size_t _i, size_t _j, size_t _k, size_t _l, size_t _value) : i(_i), j(_j), k(_k), l(_l), value(_value)
	{
		
	}

	friend std::ostream& operator << (std::ostream& of, const Cluster& cluster)
	{
		of << cluster.i << " " << cluster.j << " " << cluster.k << " " << cluster.l << " " << cluster.value << std::endl;
		return of;
	}

	bool AreClusterOverlaps(const Cluster& cluster) const
	{
		const bool x_overlap = i <= cluster.k && k >= cluster.i;
		const bool y_overlap = j <= cluster.l && l >= cluster.j;
		return x_overlap && y_overlap;
	}

	bool AreFullyOverlapped(const Cluster& cluster) const
	{
		const bool x_overlap = i <= cluster.i && k >= cluster.k;
		const bool y_overlap = j <= cluster.j && l >= cluster.l;
		return x_overlap && y_overlap;
	}

	bool operator==(const Cluster& c) const
	{
		return i == c.i && j == c.j && k == c.k && l == c.l;
	}
};

struct Hasher
{
	constexpr size_t operator()(const Cluster& c) const
	{
		return c.i * 1000 + c.j * 100 + c.k * 10 + c.l;
	}
};

constexpr std::array<Rectangle, NUMBER_OF_RECTANGLES> rectangles
{
	    Rectangle(7, 5, 10000),
		Rectangle(6, 5, 8200),
		Rectangle(5, 5, 6650),
		Rectangle(4, 5, 5250),
		Rectangle(7, 4, 5350),
		Rectangle(6, 4, 4550),
		Rectangle(5, 4, 3250),
		Rectangle(4, 4, 1850),
		Rectangle(3, 5, 1100),
		Rectangle(3, 4, 875),
		Rectangle(7, 3, 975),
		Rectangle(6, 3, 800),
		Rectangle(5, 3, 650),
		Rectangle(4, 3, 455),
		Rectangle(3, 3, 230),
		Rectangle(2, 5, 210),
		Rectangle(2, 4, 110),
		Rectangle(2, 3, 57),
		Rectangle(7, 2, 90),
		Rectangle(6, 2, 75),
		Rectangle(5, 2, 62),
		Rectangle(4, 2, 44),
		Rectangle(3, 2, 26),
		Rectangle(2, 2, 11),
		Rectangle(1, 5, 9),
		Rectangle(1, 4, 7),
		Rectangle(1, 3, 5),
		Rectangle(1, 2, 3),
		Rectangle(7, 1, 6),
		Rectangle(6, 1, 5),
		Rectangle(5, 1, 4),
		Rectangle(4, 1, 3),
		Rectangle(3, 1, 2),
		Rectangle(2, 1, 1),
		Rectangle(1, 1, 0),
};

size_t GetRandomNumber(size_t min, size_t max)
{
	static std::random_device device{};
	static std::mt19937 engine(device());

	std::uniform_int_distribution<size_t> dist(min, max);
	return dist(engine);
}

using MATRIX = std::array<std::array<size_t, M>, N>;
MATRIX matrix;
static const std::map<size_t, Rectangle, std::greater<>> elements
{
	{ rectangles.at(0).m_value,  rectangles.at(0)  },
	{ rectangles.at(1).m_value,  rectangles.at(1) },
	{ rectangles.at(2).m_value,  rectangles.at(2) },
	{ rectangles.at(3).m_value,  rectangles.at(3) },
	{ rectangles.at(4).m_value,  rectangles.at(4) },
	{ rectangles.at(5).m_value,  rectangles.at(5) },
	{ rectangles.at(6).m_value,  rectangles.at(6) },
	{ rectangles.at(7).m_value,  rectangles.at(7) },
	{ rectangles.at(8).m_value,  rectangles.at(8) },
	{ rectangles.at(9).m_value,  rectangles.at(9) },
	{ rectangles.at(10).m_value, rectangles.at(10) },
	{ rectangles.at(11).m_value, rectangles.at(11) },
	{ rectangles.at(12).m_value, rectangles.at(12) },
	{ rectangles.at(13).m_value, rectangles.at(13) },
	{ rectangles.at(14).m_value, rectangles.at(14) },
	{ rectangles.at(15).m_value, rectangles.at(15) },
	{ rectangles.at(16).m_value, rectangles.at(16) },
	{ rectangles.at(17).m_value, rectangles.at(17) },
	{ rectangles.at(18).m_value, rectangles.at(18) },
	{ rectangles.at(19).m_value, rectangles.at(19) },
	{ rectangles.at(20).m_value, rectangles.at(20) },
	{ rectangles.at(21).m_value, rectangles.at(21) },
	{ rectangles.at(22).m_value, rectangles.at(22) },
	{ rectangles.at(23).m_value, rectangles.at(23) },
	{ rectangles.at(24).m_value, rectangles.at(24) },
	{ rectangles.at(25).m_value, rectangles.at(25) },
	{ rectangles.at(26).m_value, rectangles.at(26) },
	{ rectangles.at(27).m_value, rectangles.at(27) },
	{ rectangles.at(28).m_value, rectangles.at(28) },
	{ rectangles.at(29).m_value, rectangles.at(29) },
	{ rectangles.at(30).m_value, rectangles.at(30) },
	{ rectangles.at(31).m_value, rectangles.at(31) },
	{ rectangles.at(32).m_value, rectangles.at(32) },
	{ rectangles.at(33).m_value, rectangles.at(33) },
	{ rectangles.at(34).m_value, rectangles.at(34) },
};

std::vector<MATRIX> tests;

void PrintMatrix(const MATRIX& matrix)
{
	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < M; j++)
		{
			std::cout << matrix.at(i).at(j) << " ";
		}
		std::cout << std::endl;
	}
}

struct PairHasher
{
	bool operator() (const std::pair<size_t, size_t>& c) const
	{
		return c.first * N + c.second;
	}
};

bool comparator(const std::pair<size_t, size_t>& c, const std::pair<size_t, size_t>& c1)
{
	return true;
}


std::vector<std::unordered_set<std::pair<size_t, size_t>, PairHasher>> GetSegments(const MATRIX& matrix)
{
	std::vector<std::unordered_set<std::pair<size_t, size_t>, PairHasher>> segments;

	std::unordered_set<std::pair<size_t, size_t>, PairHasher> segment;
	std::unordered_set<std::pair<size_t, size_t>, PairHasher> found;
	std::deque<std::pair<size_t, size_t>> nodes;
	
	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < M; j++)
		{
			segment.clear();
			if (matrix[i][j] != 0 && !found.contains(std::make_pair(i, j)))
			{
				nodes.emplace_back(i, j);
				while (!nodes.empty())
				{
					auto node = nodes.front();
					if (node.first == 0 && node.second == 2)
					{
						int x = 0;
						x++;
					}
					nodes.pop_front();
					found.insert(node);
					segment.insert(std::make_pair(node.first, node.second));
					if (node.first + 1 < N && matrix[node.first + 1][node.second] == 1 && !found.contains(std::make_pair(node.first + 1, node.second)))
						nodes.emplace_back(node.first + 1, node.second);
					if (node.second + 1 < M && matrix[node.first][node.second + 1] == 1 && !found.contains(std::make_pair(node.first, node.second + 1)))
						nodes.emplace_back(node.first, node.second + 1);
					if (node.second && node.second - 1 >= 0 && matrix[node.first][node.second - 1] == 1 && !found.contains(std::make_pair(node.first, node.second - 1)))
						nodes.emplace_back(node.first, node.second - 1);
					if (node.first && node.first - 1 >= 0 && matrix[node.first - 1][node.second] == 1 && !found.contains(std::make_pair(node.first - 1, node.second)))
						nodes.emplace_back(node.first - 1, node.second);
				}
			}
			if (!segment.empty())
				segments.push_back(segment);
		}
	}

	return segments;
};


bool IsNonOverlapping(const std::vector<Cluster>& rectangles)
{
	if (rectangles.empty())
		return false;

	for (size_t i = 0; i < rectangles.size(); ++i) 
	{
		for (size_t j = i + 1; j < rectangles.size(); ++j) 
		{
			if (rectangles[i].AreClusterOverlaps(rectangles[j])) 
			{
				return false;
			}
		}
	}
	return true;
}

void generateNonOverlappingCombinations(std::vector<std::vector<Cluster>>& result,
	const std::vector<Cluster>& rectangles,
	std::vector<Cluster>& currentCombination,
	size_t startIndex, size_t& max_value)
{
	if (startIndex == rectangles.size()) 
	{
		if (IsNonOverlapping(currentCombination))
		{
			result.push_back(currentCombination);
			if (currentCombination.size() > max_value)
				max_value = currentCombination.size();
		}
		return;
	}

	// Include the current rectangle and move to the next one
	currentCombination.push_back(rectangles[startIndex]);
	generateNonOverlappingCombinations(result, rectangles, currentCombination, startIndex + 1, max_value);
	currentCombination.pop_back(); // Backtrack

	// Skip the current rectangle and move to the next one
	generateNonOverlappingCombinations(result, rectangles, currentCombination, startIndex + 1, max_value);
}
//
//size_t GetBiggestRectangles(const MATRIX& matrix, std::vector<Cluster>& found)
//{
//	std::map<size_t, std::vector<Cluster>, std::greater<>> clusters;
//	auto copy = MATRIX(matrix);
//
//	for (const auto& element : elements)
//	{
//		auto rect = element.second;
//		for (size_t i = 0; i < N; i++)
//		{
//			for (size_t j = 0; j < M; j++)
//			{
//				size_t max_i = N - i;
//				size_t max_j = M - j;
//
//				if (rect.m_area > max_i * max_j)
//					continue;
//
//				if (max_i < rect.m_N)
//					continue;
//
//				if (max_j < rect.m_M)
//					continue;
//
//
//				bool IsOkey = true;
//
//				for (size_t k = i; k < i + rect.m_N && IsOkey; k++)
//				{
//					for (size_t l = j; l < j + rect.m_M && IsOkey; l++)
//					{
//						IsOkey = IsOkey && (copy[k][l] == 1);
//					}
//				}
//
//				if (IsOkey)
//				{
//					if (clusters.contains(element.first))
//						clusters[element.first].push_back(Cluster(i, j, i + rect.m_N - 1, j + rect.m_M - 1, element.first));
//					else
//					{
//						clusters[element.first] = {};
//						clusters[element.first].push_back(Cluster(i, j, i + rect.m_N - 1, j + rect.m_M - 1, element.first));
//					}
//				}
//			}
//		}
//	}
//
//	std::map<size_t, std::vector<std::vector<Cluster>>, std::greater<>> max_potential_values;
//
//	for (const auto& cluster : clusters)
//		max_potential_values[cluster.first * cluster.second.size()].push_back(cluster.second);
//
//	std::map<size_t, std::vector<std::vector<Cluster>>, std::greater<>> possible_vectors;
//	std::vector<Cluster> additionalClusters;
//	for (const auto& specific_val : max_potential_values)
//	{
//		bool shouldBreak = false;
//		for (const auto& potential : specific_val.second)
//		{
//			std::set<Cluster, Hasher> overlapping_clusters;
//			for (size_t i = 0; i < potential.size(); i++)
//			{
//				for (size_t j = i + 1; j < potential.size(); j++)
//				{
//
//					if (potential.at(i).AreClusterOverlaps(potential.at(j)))
//					{
//						overlapping_clusters.insert(potential.at(i));
//						overlapping_clusters.insert(potential.at(j));
//					}
//				}
//			}
//			if (overlapping_clusters.empty())
//			{
//				if (!possible_vectors.empty() && possible_vectors.begin()->first > specific_val.first)
//				{
//					if (true)
//					{
//						// throw implementation if specific plus possible is greater than max
//					}
//					else
//					{
//						shouldBreak = true;
//						break;
//					}
//				}
//				for (const auto& cluster : potential)
//				{
//					found.push_back(cluster);
//					for (size_t i = cluster.i; i <= cluster.k; i++)
//					{
//						for (size_t j = cluster.j; j <= cluster.l; j++)
//						{
//							copy[i][j] = 0;
//						}
//					}
//				}
//
//				return specific_val.first + GetBiggestRectangles(copy, found);
//			}
//			else
//			{
//				size_t max_length = 0;
//				std::vector<std::vector<Cluster>> nonOverlappingCombinations;
//				std::vector<Cluster> currentCombination;
//				generateNonOverlappingCombinations(nonOverlappingCombinations, potential, currentCombination, 0, max_length);
//
//				std::vector<std::vector<Cluster>> max_values;
//				for (const auto& potential_non_overlaping : nonOverlappingCombinations)
//				{
//					if (potential_non_overlaping.size() == max_length)
//						max_values.push_back(potential_non_overlaping);
//				}
//
//				if (!max_values.empty())
//				{
//					possible_vectors[max_length * (specific_val.first / potential.size())] = max_values;
//				}
//
//			}
//		}
//		if (shouldBreak)
//		{
//			break;
//		}
//	}
//
//	for (const auto& pVec : possible_vectors)
//	{
//		int max_value_found = -1;
//		std::vector<Cluster> toReturn;
//		for (const auto& combinations : pVec.second)
//		{
//			auto _copy = MATRIX(copy);
//			for (const auto& clusterUnique : combinations)
//			{
//				for (size_t i = clusterUnique.i; i <= clusterUnique.k; i++)
//					for (size_t j = clusterUnique.j; j <= clusterUnique.l; j++)
//						_copy[i][j] = 0;
//			}
//
//			std::vector<Cluster> foundOnes;
//			int result = GetBiggestRectangles(_copy, foundOnes);
//
//			if (result > max_value_found)
//			{
//				max_value_found = result;
//				toReturn = std::vector(combinations);
//			}
//
//		}
//
//		if (max_value_found > 0)
//		{
//			for (const auto& clus : toReturn)
//			{
//				found.push_back(clus);
//				for (size_t i = clus.i; i <= clus.k; i++)
//				{
//					for (size_t j = clus.j; j <= clus.l; j++)
//					{
//						copy[i][j] = 0;
//					}
//				}
//			}
//
//			
//			return pVec.first + GetBiggestRectangles(copy, found);
//		}
//	}
//
//	return 0;
//}



size_t GetBiggestRectanglesWithSegments(const MATRIX& matrix, std::vector<Cluster>& found)
{
	auto segments = GetSegments(matrix);
	std::vector<std::map<size_t, std::vector<Cluster>, std::greater<>>> found_values;

	for (const auto& segment : segments)
	{
		std::map<size_t, std::vector<Cluster>, std::greater<>> found_values_per_segment;

		for (const auto& rectangle : rectangles)
		{
			if (rectangle.m_area > segment.size())
				continue;

			for (const auto& node : segment)
			{
				if (!segment.contains(std::make_pair(node.first + rectangle.m_N - 1, node.second)) 
					|| !segment.contains(std::make_pair(node.first, node.second + rectangle.m_M - 1)))
					continue;

				bool should_break = false;
				for (size_t i = node.first; i < node.first + rectangle.m_N && !should_break; i++)
				{
					for (size_t j = node.second; j < node.second + rectangle.m_M && !should_break; j++)
					{
						if (!segment.contains(std::make_pair(i, j)))
							should_break = true;
					}
				}

				if (!should_break)
				{
					if (found_values_per_segment.contains(rectangle.m_value))
						found_values_per_segment[rectangle.m_value].push_back(Cluster(node.first, node.second, node.first + rectangle.m_N - 1, node.second + rectangle.m_M - 1, rectangle.m_value));
					else
						found_values_per_segment[rectangle.m_value] = { Cluster(node.first, node.second, node.first + rectangle.m_N - 1, node.second + rectangle.m_M - 1, rectangle.m_value) };
				}
			}
		}

		if (!found_values_per_segment.empty())
			found_values.push_back(found_values_per_segment);
	}

	size_t result = 0;

	int clust = -1;
	for (const auto& values : found_values)
	{
		clust++;
		std::map<size_t, std::vector<Cluster>, std::greater<>> max_possible_values;
		for (const auto& max : values)
		{
			if (!max_possible_values.contains(max.first * max.second.size()))
				max_possible_values[max.first * max.second.size()] = std::vector<Cluster>();

			for (const auto& cluster : max.second)
			{
				max_possible_values[max.first * max.second.size()].push_back(cluster);
			}
		}

		if (!max_possible_values.empty() && (*max_possible_values.begin()).first == 0)
		{

			const auto& clusters = (*max_possible_values.begin()).second;

			for (const auto& clusss : clusters)
				found.push_back(clusss);

			continue;
		}

		std::map < size_t, std::vector<std::vector<Cluster>>, std::greater<>> biggest_combinations;

		for (const auto& try_val : max_possible_values)
		{
			if (try_val.second.size() == 1)
			{
				if (!biggest_combinations.empty() && (*biggest_combinations.begin()).first > try_val.first)
					break;

				const auto& clus = try_val.second.front();
				found.push_back(try_val.second.front());
				MATRIX _matrix;
				for (size_t i = 0; i < N; i++)
				{
					for (size_t j = 0; j < M; j++)
					{
						if (i >= clus.i && i <= clus.k && j >= clus.j && j <= clus.l)
							_matrix[i][j] = 0;
						else if (segments[clust].contains(std::make_pair(i, j)))
							_matrix[i][j] = 1;
						else
							_matrix[i][j] = 0;
					}
				}

				result += (try_val.first + GetBiggestRectanglesWithSegments(_matrix, found));
				biggest_combinations.clear();
				break;
			}
			else
			{
				if (!biggest_combinations.empty() && (*biggest_combinations.begin()).first > try_val.first)
					break;
				const auto& elements_to_pick = try_val.second;
				size_t max_length = 0;
				std::vector<std::vector<Cluster>> nonOverlappingCombinations;
				std::vector<Cluster> currentCombination;
				generateNonOverlappingCombinations(nonOverlappingCombinations, elements_to_pick, currentCombination, 0, max_length);

				std::vector<std::vector<Cluster>> max_values;
				for (const auto& potential_non_overlaping : nonOverlappingCombinations)
				{
					if (potential_non_overlaping.size() == max_length)
						max_values.push_back(potential_non_overlaping);
				}

				if (!max_values.empty())
				{
					biggest_combinations[max_length * (try_val.first / elements_to_pick.size())] = max_values;
				}
			}
		}

		if (!biggest_combinations.empty())
		{
			const auto& values_to_check = *biggest_combinations.begin();
			
			int max = -1;
			std::vector<Cluster> biggest;

			for (const auto& specific : values_to_check.second)
			{
				MATRIX _mmatrix;
				for (const auto& spec : specific)
				{
					for (size_t i = 0; i < N; i++)
					{
						for (size_t j = 0; j < M; j++)
						{
							if (i >= spec.i && i <= spec.k && j >= spec.j && j <= spec.l)
								_mmatrix[i][j] = 0;
							else if (segments[clust].contains(std::make_pair(i, j)))
								_mmatrix[i][j] = 1;
							else
								_mmatrix[i][j] = 0;
						}
					}
				}

				std::vector<Cluster> foundsss;
				auto resultss = GetBiggestRectanglesWithSegments(_mmatrix, foundsss);

				if (int(resultss) > max)
				{
					max = resultss;
					biggest = std::vector(specific);
				}

			}

			if (!biggest.empty())
			{
				MATRIX _mmatrix;
				
				for (size_t i = 0; i < N; i++)
				{
					for (size_t j = 0; j < M; j++)
					{
						if (segments[clust].contains(std::make_pair(i, j)))
							_mmatrix[i][j] = 1;
						else
							_mmatrix[i][j] = 0;
					}
				}

				for (size_t i = 0; i < N; i++)
				{
					for (size_t j = 0; j < M; j++)
					{
						for (const auto& spec : biggest)
						{
							if (i >= spec.i && i <= spec.k && j >= spec.j && j <= spec.l)
								_mmatrix[i][j] = 0;
						}
					}
				}

				for (const auto& spec : biggest)
				{	
					found.push_back(spec);
				}
				result += (values_to_check.first + GetBiggestRectanglesWithSegments(_mmatrix, found));
			}

		}

	}

	return result;
}


size_t GetBiggestRectanglesEasy(MATRIX& matrix, std::vector<Cluster>& found)
{
	for (const auto& element : rectangles)
	{
		std::vector<Cluster> clusters;

		const auto& rect = element;
		for (size_t i = 0; i < N; i++)
		{
			for (size_t j = 0; j < M; j++)
			{
				size_t max_i = N - i;
				size_t max_j = M - j;

				if (rect.m_area > max_i * max_j)
					continue;

				if (max_i < rect.m_N)
					continue;

				if (max_j < rect.m_M)
					continue;


				bool IsOkey = true;

				for (size_t k = i; k < i + rect.m_N && IsOkey; k++)
				{
					for (size_t l = j; l < j + rect.m_M && IsOkey; l++)
					{
						IsOkey = IsOkey && (matrix[k][l] == 1);
					}
				}

				if (IsOkey)
				{
					clusters.push_back(Cluster(i, j, i + rect.m_N - 1, j + rect.m_M - 1, element.m_value));
				}
			}
		}

		if (clusters.size() == 1)
		{
			const auto& cluster = clusters.front();
			found.push_back(cluster);
			for (size_t i = cluster.i; i <= cluster.k; i++)
			{
				for (size_t j = cluster.j; j <= cluster.l; j++)
				{
					matrix[i][j] = 0;
				}
			}
			return element.m_value + GetBiggestRectanglesEasy(matrix, found);
		}
		else
		{
			const auto& elements_to_pick = clusters;
			size_t max_length = 0;
			std::vector<std::vector<Cluster>> nonOverlappingCombinations;
			std::vector<Cluster> currentCombination;
			generateNonOverlappingCombinations(nonOverlappingCombinations, elements_to_pick, currentCombination, 0, max_length);

			std::vector<std::vector<Cluster>> max_values;
			for (const auto& potential_non_overlaping : nonOverlappingCombinations)
			{
				if (potential_non_overlaping.size() == max_length)
					max_values.push_back(potential_non_overlaping);
			}

			int max = -1;
			std::vector<Cluster> strongest;
			for (const auto& values : max_values)
			{
				MATRIX copy = MATRIX(matrix);
				std::vector<Cluster> foundCopy;

				for (const auto& small_rect : values)
				{
					for (size_t i = small_rect.i; i <= small_rect.k; i++)
					{
						for (size_t j = small_rect.j; j <= small_rect.l; j++)
						{
							copy[i][j] = 0;
						}
					}
				}
				int result_temp = GetBiggestRectanglesEasy(copy, foundCopy);
				if (result_temp > max)
				{
					max = result_temp;
					strongest = std::vector(values);
				}
			}

			if (!strongest.empty())
			{
				for (const auto& small_rect : strongest)
				{
					for (size_t i = small_rect.i; i <= small_rect.k; i++)
					{
						for (size_t j = small_rect.j; j <= small_rect.l; j++)
						{
							matrix[i][j] = 0;
						}
					}
					found.push_back(small_rect);
				}
				return element.m_value * strongest.size() + GetBiggestRectanglesEasy(matrix, found);
			}
		}
	}
	return 0;
}


void CreateMatrix(long code, MATRIX& x)
{
	int shift = 0;
	for (int j = 0; j < 7; j++)
	{
		for (int i = 0; i < 5; i++)
		{
			if ((i > 0 && i < 4) && (j == 2 || j == 4))
			{
				x[i][j] = 1;
			}
			else
			{
				x[i][j] = ((code >> shift) & 1);
				shift++;
			}
		}
	}
}


int main()
{
	std::ifstream file("Tests.txt");
	struct TestValues
	{
		MATRIX matrix;
		size_t value;
		size_t code;
	};
	unsigned long long tests = 0;
	double time = 0;
	for (auto path : std::filesystem::directory_iterator(R"(C:\Users\radoi\Desktop\wintables_corrected)"))
	{
		auto fileName = path.path().filename().string();
		auto winType = fileName.find('_');
		std::string copy = std::string(fileName);
		copy.erase(0, winType + 1);
		std::string number = "";
		auto c = 0;
		while (std::isdigit(copy[c]) && c < copy.size())
		{
			number.push_back(copy[c]);
			c++;
		}

		std::ifstream open(path.path().string().c_str());
		std::vector<long> values;
		std::vector<TestValues> test_values;
		if (open.is_open())
		{
			std::string line;
			while (std::getline(open, line))
			{
				std::string number;
				int last = line.size() - 1;
				while (last >= 0 && std::isdigit(line[last]))
				{
					number.push_back(line[last]);
					last--;
				}

				std::reverse(number.begin(), number.end());
				values.push_back(atol(number.c_str()));
			}
			open.close();
		}
		for (const auto& val : values)
		{
			TestValues val_;
			MATRIX m;
			CreateMatrix(val, m);
			val_.code = val;
			val_.matrix = MATRIX(m);
			val_.value = atol(number.c_str());
			test_values.push_back(val_);
		}
		//std::cout << "Begin test for : " << atol(number.c_str()) << std::endl;

		auto start = std::chrono::high_resolution_clock::now();

		size_t failed = 0;
		size_t success = 0;

		for (auto& test : test_values)
		{
			tests++;
			MATRIX copy = MATRIX(test.matrix);
			std::vector<Cluster> clusters_found;
			if (auto result = GetBiggestRectanglesEasy(copy, clusters_found); result != test.value)
			{
				failed++;
				std::cout << "Wrong results" << std::endl;
				PrintMatrix(test.matrix);
				std::cout << "Expected: " << test.value << " , got: " << result << std::endl;
				for (const auto& cluster : clusters_found)
					std::cout << cluster << std::endl;
				clusters_found.clear();
				copy = MATRIX(test.matrix);
				auto cc = GetBiggestRectanglesEasy(copy, clusters_found);

			}
			else
			{
				success++;
			}
		}
		auto end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> duration = (end - start);
		time += duration.count();
		//std::cout << "Number of lines" << values.size() << std::endl;
		//std::cout << "End test for : " << atol(number.c_str()) << std::endl;
	}

	std::cout << "Passed tests: " << tests << std::endl;
	std::cout << "Time: " << time << std::endl;

}