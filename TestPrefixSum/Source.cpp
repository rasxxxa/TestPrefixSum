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

void GenerateNonOverlappingCombinations(std::vector<std::vector<Cluster>>& result,
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
	GenerateNonOverlappingCombinations(result, rectangles, currentCombination, startIndex + 1, max_value);
	currentCombination.pop_back(); // Backtrack

	// Skip the current rectangle and move to the next one
	GenerateNonOverlappingCombinations(result, rectangles, currentCombination, startIndex + 1, max_value);
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
			GenerateNonOverlappingCombinations(nonOverlappingCombinations, elements_to_pick, currentCombination, 0, max_length);

			std::vector<std::vector<Cluster>> max_values;
			for (const auto& potential_non_overlaping : nonOverlappingCombinations)
			{
				if (potential_non_overlaping.size() == max_length)
					max_values.push_back(potential_non_overlaping);
			}

			int max = -1;
			std::vector<Cluster>* strongest = nullptr;
			for (auto& values : max_values)
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
					strongest = (&values);
				}
			}

			if (strongest && !(*strongest).empty())
			{
				for (const auto& small_rect : (*strongest))
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
				return element.m_value * (*strongest).size() + GetBiggestRectanglesEasy(matrix, found);
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
	////std::ifstream file("Tests.txt");
	struct TestValues
	{
		MATRIX matrix;
		size_t value;
		size_t code;
	};
	unsigned long long tests = 0;
	double time = 0;
	for (auto path : std::filesystem::directory_iterator(R"(C:\wintables_corrected)"))
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
		auto start = std::chrono::high_resolution_clock::now();

		size_t failed = 0;
		size_t success = 0;

		for (auto& test : test_values)
		{
			tests++;
			std::vector<Cluster> clusters_found;
			if (auto result = GetBiggestRectanglesEasy(test.matrix, clusters_found); result != test.value)
			{
				failed++;
				std::cout << "Wrong results" << std::endl;
				PrintMatrix(test.matrix);
				std::cout << "Expected: " << test.value << " , got: " << result << std::endl;
				for (const auto& cluster : clusters_found)
					std::cout << cluster << std::endl;
				clusters_found.clear();

			}
			else
			{
				success++;
			}
		}
		auto end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> duration = (end - start);
		time += duration.count();
		break;
	}

	std::cout << "Passed tests: " << tests << std::endl;
	std::cout << "Time: " << time << std::endl;
}