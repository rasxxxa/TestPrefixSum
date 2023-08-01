#include <iostream>
#include <array>
#include <random>
#include <fstream>
#include <map>
#include <functional>
#include <vector>
#include <chrono>
#include <unordered_set>

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

size_t GetBiggestRectangles(const MATRIX& matrix, std::vector<Cluster>& found)
{
	std::map<size_t, std::vector<Cluster>, std::greater<>> clusters;
	MATRIX copy = MATRIX(matrix);

	size_t iteration = 0;
	size_t returnValue = 0;

	for (const auto& element : elements)
	{
		auto rect = element.second;
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
						IsOkey = IsOkey && (copy[k][l] == 1);
						iteration++;
					}
				}

				if (IsOkey)
				{
					if (clusters.contains(element.first))
						clusters[element.first].push_back(Cluster(i, j, i + rect.m_N - 1, j + rect.m_M - 1, element.first));
					else
						clusters[element.first] = { Cluster(i, j, i + rect.m_N - 1, j + rect.m_M - 1, element.first) };
				}
			}
		}
	}

	for (auto it = clusters.begin(); it != clusters.end(); ++it)
	{
		if (clusters.contains(it->first))
		{
			const auto cluster = clusters[it->first];
			if (cluster.size() == 1)
			{
				const auto& clusterUnique = cluster.front();
				for (size_t i = clusterUnique.i; i <= clusterUnique.k; i++)
					for (size_t j = clusterUnique.j; j <= clusterUnique.l; j++)
						copy[i][j] = 0;

				returnValue += clusterUnique.value;
				found.push_back(clusterUnique);
			}
			else if (!cluster.empty())
			{
				std::unordered_set<Cluster, Hasher> overlapping_clusters;
				const auto& cluster_vec = cluster;
				for (size_t i = 0; i < cluster_vec.size(); i++)
				{
					for (size_t j = i; j < cluster_vec.size(); j++)
					{
						if (cluster_vec.at(i).AreClusterOverlaps(cluster_vec.at(j)))
						{
							overlapping_clusters.insert(cluster_vec.at(i));
							overlapping_clusters.insert(cluster_vec.at(j));
						}
					}
				}

				for (size_t i = 0; i < cluster_vec.size(); i++)
				{
					if (overlapping_clusters.contains(cluster_vec[i]))
					{
						const auto& clusterUnique = cluster_vec[i];
						for (size_t i = clusterUnique.i; i < clusterUnique.k; i++)
							for (size_t j = clusterUnique.j; j < clusterUnique.l; j++)
								copy[i][j] = 0;
						returnValue += clusterUnique.value;
					}
				}

				if (!overlapping_clusters.empty())
				{
					std::vector<Cluster> _rectangles;
					for (const auto& rectangle : overlapping_clusters)
						_rectangles.push_back(rectangle);

					size_t max_length = 0;
					std::vector<std::vector<Cluster>> nonOverlappingCombinations;
					std::vector<Cluster> currentCombination;

					generateNonOverlappingCombinations(nonOverlappingCombinations, _rectangles, currentCombination, 0, max_length);

					int max_value_found = INT_MIN;
					int index = -1;

					for (size_t rectangle = 0; rectangle < nonOverlappingCombinations.size(); rectangle++)
					{
						if (nonOverlappingCombinations[rectangle].size() != max_length)
							continue;

						auto newCopy = MATRIX(copy);
						for (const auto& clusterToRemove : nonOverlappingCombinations[rectangle])
						{
							for (size_t ii = clusterToRemove.i; ii <= clusterToRemove.k; ii++)
								for (size_t jj = clusterToRemove.j; jj <= clusterToRemove.l; jj++)
									newCopy[ii][jj] = 0;

						}
						std::vector<Cluster> clustersToRemove;
						size_t result = GetBiggestRectangles(newCopy, clustersToRemove);
							
						if (const int castedValue = static_cast<int>(result); castedValue > max_value_found && result > 0)
						{
							max_value_found = castedValue;
							index = static_cast<int>(rectangle);
						}

					}

					if (max_value_found == INT_MIN)
					{
						// ASSERT error;

					}
					else
					{
						for (const auto& clusterToRemove : nonOverlappingCombinations[index])
						{
							for (size_t ii = clusterToRemove.i; ii <= clusterToRemove.k; ii++)
								for (size_t jj = clusterToRemove.j; jj <= clusterToRemove.l; jj++)
									copy[ii][jj] = 0;

							returnValue += clusterToRemove.value;
							found.push_back(clusterToRemove);
						}
					}

				}
			}
		}
	}
	return returnValue;
}




int main()
{
	MATRIX matrix;
	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < M; j++)
		{
			matrix[i][j] = 1;
			if (i == 0)
				matrix[i][j] = 0;
		}
	}
	matrix[4][3] = 0;
	std::cout << std::endl;
	//matrix[0][0] = 1;
	//matrix[0][1] = 1;
	//matrix[0][3] = 1;
	//matrix[0][4] = 1;
	//matrix[0][5] = 1;
	//matrix[1][0] = 1;
	//matrix[1][1] = 1;
	//matrix[1][2] = 1;
	//matrix[1][4] = 1;
	//matrix[1][5] = 1;
	//matrix[1][6] = 1;
	//matrix[2][0] = 1;
	//matrix[2][1] = 1;
	//matrix[2][2] = 1;
	//matrix[2][3] = 1;
	//matrix[2][5] = 1;
	//matrix[2][6] = 1;
	//matrix[3][1] = 1;
	//matrix[3][2] = 1;
	//matrix[3][3] = 1;
	//matrix[4][2] = 1;
	//matrix[4][3] = 1;

	PrintMatrix(matrix);
	std::vector<Cluster> clusters_found;
	auto start = std::chrono::high_resolution_clock::now();
	std::cout << "Value " << GetBiggestRectangles(matrix, clusters_found) << std::endl;
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> duration = end - start;


	std::cout << "Duration needed: " << duration.count() << " seconds" << std::endl;
	for (const auto& cluster : clusters_found)
	{
		std::cout << "Value: " << cluster.value << std::endl;
		std::cout << cluster << std::endl;
	}
}