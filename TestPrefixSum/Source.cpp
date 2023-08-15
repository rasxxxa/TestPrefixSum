#include <iostream>
#include <array>
#include <random>
#include <fstream>
#include <map>
#include <functional>
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

	for (const auto& values : found_values)
	{
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

		std::map < size_t, std::vector<std::vector<Cluster>>, std::greater<>> biggest_combinations;

		for (const auto& try_val : max_possible_values)
		{
			if (try_val.second.size() == 1)
			{
				const auto& clus = try_val.second.front();
				found.push_back(try_val.second.front());
				MATRIX _matrix;
				for (size_t i = 0; i < N; i++)
				{
					for (size_t j = 0; j < M; j++)
					{
						if (i >= clus.i && i <= clus.k && j >= clus.j && j <= clus.l)
							_matrix[i][j] = 1;
						else
							_matrix[i][j] = 0;
					}
				}

				result += (try_val.first + GetBiggestRectanglesWithSegments(_matrix, found));
				break;
			}
			else
			{
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
			// todo
		}

	}

	return result;
}



int main()
{
	std::ifstream file("Tests.txt");
	struct TestValues
	{
		MATRIX matrix;
		size_t value;
	};


	std::vector<TestValues> test_values;
	if (file.is_open())
	{
		std::string line;
		std::vector<std::string> lines;

		while (std::getline(file, line))
		{
			lines.push_back(line);
			std::vector<std::vector<int>> matrix;
			if (lines.size() == N + 1)
			{
				for (int i = 0; i < N; i++)
				{
					std::vector<int> row(M, 0);
					int elements_put_in_row = 0;
					std::string specific_row = lines[i];
					for (int j = 0; j < specific_row.size(); j++)
					{
						if (specific_row[j] == '1')
						{
							row[elements_put_in_row] = 1;
							elements_put_in_row++;
						}
						else if (specific_row[j] == '0')
						{
							row[elements_put_in_row] = 0;
							elements_put_in_row++;
						}
					}
					matrix.push_back(row);
				}
				
				size_t value = atol(lines[5].c_str());

				TestValues test;
				test.value = value;

				for (int i = 0; i < N; i++)
				{
					for (int j = 0; j < M; j++)
					{
						test.matrix[i][j] = matrix[i][j];
					}
				}
				test_values.push_back(test);
				lines.clear();
			}
		}


		file.close();
	}
	
	auto start = std::chrono::high_resolution_clock::now();
	
	std::cout << "Begin test" << std::endl;
	size_t failed = 0;
	for (const auto& test : test_values)
	{
		//auto segments = GetSegments(test.matrix);

		//for (const auto& s : segments)
		//{
		//	std::cout << "segment" << std::endl;
		//	for (const auto& v : s)
		//		std::cout << "[" << v.first << "," << v.second << "]" << std::endl;
		//}

		PrintMatrix(test.matrix);
		std::vector<Cluster> clusters_found;
		auto result = GetBiggestRectanglesWithSegments(test.matrix, clusters_found);
		std::cout << "Wrong results" << std::endl;
		//PrintMatrix(test.matrix);
		std::cout << "Expected: " << test.value << " , got: " << result << std::endl;
		for (const auto& cluster : clusters_found)
			std::cout << cluster << std::endl;
		clusters_found.clear();
		//{
		//	
		//	//GetBiggestRectangles(test.matrix, clusters_found);
		//}
		//else
		//{
		//	//std::cout << "Test passed " << std::endl;
		//}
	}
	std::cout << "Failed for " << failed << std::endl; 
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> duration = end - start;

	std::cout << "Time passed for " << test_values.size() << " tests : " << duration.count() << " seconds" << std::endl;
}