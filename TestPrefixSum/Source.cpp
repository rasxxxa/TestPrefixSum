#include <iostream>
#include <array>
#include <random>
#include <functional>
struct Cluster
{
	size_t horizontal{ 0 };
	size_t vertical{ 0 };
};

constexpr size_t M = 6;
constexpr size_t N = 6;

std::array<std::array<size_t, M>, N> _matrix;
std::array<std::array<Cluster, M>, N> _prefix;

size_t GetRandomNumber(size_t min, size_t max)
{
	static std::random_device device{};
	static std::mt19937 engine(device());

	std::uniform_int_distribution<size_t> dist(min, max);
	return dist(engine);
}

void PrintMatrix(const std::array<std::array<size_t, M>, N>& matrix)
{
	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < M; j++)
		{
			std::cout << matrix[i][j] << " ";
		}
		std::cout << std::endl;
	}
}

void PrintMatrix(const std::array<std::array<Cluster, M>, N>& matrix)
{
	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < M; j++)
		{
			std::cout << "(" << matrix[i][j].horizontal  << "," << matrix[i][j].vertical<< ")" << " ";
		}
		std::cout << std::endl;
	}
}

void FillWithRandoms(std::array<std::array<size_t, M>, N>& matrix)
{
	for (size_t i = 0; i < N; i++)
		for (size_t j = 0; j < M; j++)
			matrix[i][j] = GetRandomNumber(0, 1);
}

void PrintMax(const std::array<std::array<Cluster, M>, N>& prefix)
{
	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < M; j++)
		{
			std::cout << "((" << (prefix[i][j].horizontal + 1) * (prefix[i][j].vertical + 1) << "))";
		}
		std::cout << std::endl;
	}
}

void CalculatePrefix(const std::array<std::array<size_t, M>, N>& matrix, std::array<std::array<Cluster, M>, N>& prefix)
{
	
	for (size_t i = 0; i < M; i++)
	{
		for (size_t j = 0; j < N; j++)
		{
			if (!matrix[i][j])
				continue;

			int posI = i - 1; 
			int posJ = j - 1;
			Cluster cluster{};

			if (posI >= 0)
			{
				if (matrix[posI][j])
					cluster.vertical += prefix[posI][j].vertical + 1;
			}
			if (posJ >= 0)
			{
				if (matrix[i][posJ])
					cluster.horizontal += prefix[i][posJ].horizontal + 1;
			}
			prefix[i][j] = cluster;
		}
	}
}

size_t hash_matrix(const std::array<std::array<size_t, N>, M>& matrix)
{
	size_t hash = 10;
	for (size_t i = 0; i < N; i++)
		for (size_t j = 0; j < M; j++)
			hash ^= std::hash<bool>{}(matrix[i][j]);

	return hash;
}

int main()
{
	FillWithRandoms(_matrix);
	PrintMatrix(_matrix);
	std::cout << hash_matrix(_matrix);
	/*CalculatePrefix(_matrix, _prefix);
	std::cout << std::endl << std::endl << std::endl;
	PrintMatrix(_prefix);
	std::cout << std::endl << std::endl << std::endl;
	PrintMax(_prefix);*/



}