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

	constexpr bool operator == (const Rectangle& rect) const
	{
		return m_M == rect.m_M && m_N == rect.m_N;
	}

	bool FullyContainsRectangle(const Rectangle& rect) const
	{
		return m_M >= rect.m_M && m_N >= rect.m_N;
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

	currentCombination.push_back(rectangles[startIndex]);
	GenerateNonOverlappingCombinations(result, rectangles, currentCombination, startIndex + 1, max_value);
	currentCombination.pop_back();

	GenerateNonOverlappingCombinations(result, rectangles, currentCombination, startIndex + 1, max_value);
}

void GetClusterForSpecificRect(std::vector<Cluster>& clusters, const MATRIX& matrix, const Rectangle& rect)
{
	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < M; j++)
		{

			size_t max_i = N - i;
			size_t max_j = M - j;

			if (!((rect.m_area > max_i * max_j) || (max_i < rect.m_N) || (max_j < rect.m_M)))
			{
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
					clusters.emplace_back(i, j, i + rect.m_N - 1, j + rect.m_M - 1, rect.m_value);
				}
			}
		}
	}
}

void InitializeMatrixWithSpecificRect(MATRIX& matrix, const Cluster& cluster, size_t value)
{
	for (size_t i = cluster.i; i <= cluster.k; i++)
		for (size_t j = cluster.j; j <= cluster.l; j++)
			matrix[i][j] = value;
}

size_t GetBiggestRectanglesEasy(size_t initialDepth, MATRIX& matrix, std::vector<Cluster>& found)
{
	const auto& rect = rectangles[initialDepth];
	if (rect.m_area == 1)
	{
		for (size_t row = 0; row < N; row++)
			for (size_t column = 0; column < M; column++)
				if (matrix[row][column])
					found.emplace_back(row, column, row, column, rect.m_value);

		return 0;
	}

	std::vector<Cluster> clusters;
	GetClusterForSpecificRect(clusters, matrix, rect);

	if (clusters.empty())
	{
		return GetBiggestRectanglesEasy(initialDepth + 1, matrix, found);
	}
	else if (clusters.size() == 1)
	{
		const auto& cluster = clusters.front();
		found.push_back(cluster);
		InitializeMatrixWithSpecificRect(matrix, cluster, 0);
		return rect.m_value + GetBiggestRectanglesEasy(initialDepth + 1, matrix, found);
	}
	else
	{
		const auto& elements_to_pick = clusters;
		size_t max_length = 0;
		std::vector<std::vector<Cluster>> nonOverlappingCombinations;
		std::vector<Cluster> currentCombination;
		GenerateNonOverlappingCombinations(nonOverlappingCombinations, elements_to_pick, currentCombination, 0, max_length);

		int max = -1;
		std::vector<Cluster>* strongest = nullptr;
		std::vector<Cluster> foundInner;
		for (auto& values : nonOverlappingCombinations)
		{
			if (values.size() == max_length)
			{
				MATRIX copy = MATRIX(matrix);
				std::vector<Cluster> foundCopy;

				for (const auto& small_rect : values)
					InitializeMatrixWithSpecificRect(copy, small_rect, 0);

				int result_temp = GetBiggestRectanglesEasy(initialDepth + 1, copy, foundCopy);
				if (result_temp > max)
				{
					max = result_temp;
					strongest = (&values);
					foundInner = foundCopy;
				}
			}
		}

		if (strongest && !(*strongest).empty())
		{
			found.insert(found.end(), (*strongest).begin(), (*strongest).end());
			found.insert(found.end(), foundInner.begin(), foundInner.end());
			return rect.m_value * (*strongest).size() + max;
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


//consteval std::array<Rectangle, NUMBER_OF_RECTANGLES> RectangleToCheck(const Rectangle& minimal, size_t& returned)
//{
//	constexpr std::array<Rectangle, NUMBER_OF_RECTANGLES> to_check(rectangles);
//	auto first = std::ranges::find(rectangles, minimal);
//	returned = first->m_area;
//	return to_check;
//}

std::vector<Rectangle> RectanglesToCheck(const Rectangle& minimal)
{
	std::vector<Rectangle> to_check;
	std::vector<Rectangle> candidates;
	auto first = std::ranges::find(rectangles, minimal);
	if (first != rectangles.end())
	{
		for (auto it = rectangles.begin(); it != first; it++)
			candidates.push_back(*it);

		to_check.push_back(minimal);
		while (!candidates.empty())
		{
			auto rect_minimal = to_check.back();
			
				auto [first, last] = std::ranges::remove_if(candidates, [rect = rect_minimal](const Rectangle& rectangle) {
				
				return rectangle.FullyContainsRectangle(rect);
				
				});

			candidates.erase(first, last);
			if (!candidates.empty())
			{
				to_check.push_back(candidates.back());
				candidates.pop_back();
			}
		}

	}
	return to_check;
}

std::vector<std::pair<size_t, size_t>> PossibleFields(MATRIX& matrix)
{
	constexpr Rectangle smallest(4, 2, 1000);
	std::vector<std::pair<size_t, size_t>> teaser_fields;
	static auto wins = RectanglesToCheck(smallest);
	static constexpr auto CheckSubField = [](const MATRIX& m, const Rectangle& rect, size_t left, size_t right, size_t top, size_t bottom) 
	{
		for (int i = top; i < bottom; i++)
		{
			for (int j = left; j < right; j++)
			{
				if (m[i][j] == 0)
					continue;

				size_t max_i = N - i;
				size_t max_j = M - j;

				if (!((rect.m_area > max_i * max_j) || (max_i < rect.m_N) || (max_j < rect.m_M)))
				{
					bool found = true;

					for (int k = i; k < i + rect.m_N && found; k++)
					{
						for (int l = j; l < j + rect.m_M && found; l++)
							found = found && (m[k][l] == 1);
					}
					if (found)
						return true;

				}
			}
		}

		return false;
	};

	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < M; j++)
		{
			if (matrix[i][j])
				continue;

			matrix[i][j] = 1;
			for (const auto& rect : wins)
			{
				int top = std::max(0, int(i - rect.m_N + 1));
				int bottom = std::min(int(N), int(i + rect.m_N - 1));
				int left = std::max(0, int(j - rect.m_M + 1));
				int right = std::min(int(M), int(j + rect.m_M - 1));
				if (CheckSubField(matrix, rect, left, right, top, bottom))
				{
					teaser_fields.push_back({ i, j });
					CheckSubField(matrix, rect, left, right, top, bottom);
					break;
				}
			}

			matrix[i][j] = 0;
		}
	}

	return teaser_fields;
}


void DoTests()
{
	////std::ifstream file("Tests.txt");
	struct TestValues
	{
		MATRIX matrix;
		size_t value;
		size_t code;
	};
	double fullTime = 0.0;
	unsigned long long tests = 0;
	double maxTime = -100.0;
	long valueMax = 0;
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
		

		size_t failed = 0;
		size_t success = 0;

		auto fullStart = std::chrono::high_resolution_clock::now();
		for (auto& test : test_values)
		{
			tests++;
			std::vector<Cluster> clusters_found;
			MATRIX copy(test.matrix);
			auto start = std::chrono::high_resolution_clock::now();
			if (auto result = GetBiggestRectanglesEasy(0, test.matrix, clusters_found); result != test.value)
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
				auto res = PossibleFields(copy);
				success++;
			}
			auto end = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> duration1 = (end - start);
			if (auto dur = duration1.count(); dur > maxTime)
			{
				maxTime = dur;
				valueMax = test.code;
			}
		}
		auto fullEnd = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> duration = (fullEnd - fullStart);
		fullTime += duration.count();
		break;
	}

	std::cout << "Passed tests: " << tests << std::endl;
	std::cout << "max time passed: " << maxTime << std::endl;
	std::cout << "Most complex combination " << valueMax << std::endl;
	std::cout << "Full time: " << fullTime << std::endl;
}



int main()
{
	DoTests();
	//MATRIX m;
	//CreateMatrix(254180095, m);
	//std::vector<Cluster> clusters_found;
	//auto start = std::chrono::high_resolution_clock::now();
	//PrintMatrix(m);
	//std::cout << GetBiggestRectanglesEasy(0, m, clusters_found) << std::endl;
	//for (const auto& rect : clusters_found)
	//	std::cout << rect << std::endl;

	//auto end = std::chrono::high_resolution_clock::now();
	//std::chrono::duration<double> duration1 = (end - start);
	//std::cout << duration1.count() << std::endl;
}