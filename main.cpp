#include <iostream>
#include <iomanip>
#include <vector>
#include <numeric>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <random>

#include "chi_squar.hpp"
#include "nist_test.hpp"
#include "modified_generators.hpp"
#include "standart_generators.hpp"

/** @brief Количество выборок для статистического анализа */
constexpr size_t NUM_SAMPLES = 20;

/** @brief Размер каждой выборки (количество чисел) */
constexpr size_t SAMPLE_SIZE = 20000;

/** @brief Уровень значимости для статистических тестов (α = 0.05) */
constexpr double SIGNIFICANCE = 0.05;

/**
 * @brief Вычисляет выборочное среднее арифметическое.
 * @tparam T Тип элементов выборки
 * @param sample Вектор значений
 * @return Среднее арифметическое
 */
template <typename T>
double mean(const std::vector<T>& sample) {
    return std::accumulate(sample.begin(), sample.end(), 0.0) / sample.size();
}

/**
 * @brief Вычисляет выборочное стандартное отклонение.
 * @tparam T Тип элементов выборки
 * @param sample Вектор значений
 * @return Стандартное отклонение
 */
template <typename T>
double stddev(const std::vector<T>& sample) {
    double m = mean(sample);
    double sq = 0.0;
    for (const auto& x : sample) {
        sq += (x - m) * (x - m);
    }
    return std::sqrt(sq / sample.size());
}

/**
 * @brief Вычисляет коэффициент вариации (CV = stddev / mean).
 * @tparam T Тип элементов выборки
 * @param sample Вектор значений
 * @return Коэффициент вариации
 */
template <typename T>
double variability(const std::vector<T>& sample) {
    return stddev(sample) / mean(sample);
}

/**
 * @brief Генерирует заданное количество выборок указанного размера.
 * @tparam Generator Тип генератора (должен иметь метод next())
 * @param gen Генератор (передаётся по ссылке)
 * @param numSamples Количество выборок
 * @param sampleSize Размер каждой выборки
 * @return Вектор векторов сгенерированных чисел
 */
template <typename Generator>
std::vector<std::vector<typename Generator::result_type>>
generateSamples(Generator& gen, size_t numSamples, size_t sampleSize) {
    std::vector<std::vector<typename Generator::result_type>> samples(numSamples);
    for (size_t i = 0; i < numSamples; ++i) {
        samples[i].reserve(sampleSize);
        for (size_t j = 0; j < sampleSize; ++j) {
            samples[i].push_back(gen.next());
        }
    }
    return samples;
}

/**
 * @brief Выводит статистики и результат хи-квадрат для каждой выборки.
 * @tparam T Тип элементов выборки
 * @param name Имя генератора
 * @param samples Вектор выборок
 * @param minVal Минимальное значение генератора
 * @param maxVal Максимальное значение генератора
 * @param significance Уровень значимости
 */
template <typename T>
void printStats(const std::string& name,
                const std::vector<std::vector<T>>& samples,
                T minVal, T maxVal, double significance) {
    std::cout << "\n=== " << name << " ===\n";
    size_t passed = 0;
    for (size_t i = 0; i < samples.size(); ++i) {
        double m = mean(samples[i]);
        double sd = stddev(samples[i]);
        double cv = variability(samples[i]);
        bool chi_ok = custom::chi_squared_test(samples[i], minVal, maxVal,
                                               significance, custom::DistributionType::UNIFORM_INT);
        if (chi_ok) ++passed;
        std::cout << "Sample " << std::setw(2) << i + 1
                  << ": mean=" << std::fixed << std::setprecision(2) << m
                  << ", std=" << sd
                  << ", CV=" << cv
                  << ", Chi-square=" << (chi_ok ? "PASS" : "FAIL") << "\n";
    }
    std::cout << "Chi-square passed: " << passed << "/" << samples.size() << "\n";
}

/**
 * @brief Преобразует вектор целых чисел в битовую последовательность.
 * @tparam T Целочисленный тип
 * @param sample Вектор чисел
 * @param bitsPerNumber Сколько младших битов использовать
 * @return Вектор битов (bool)
 */
template <typename T>
std::vector<bool> toBits(const std::vector<T>& sample, int bitsPerNumber) {
    std::vector<bool> bits;
    bits.reserve(sample.size() * bitsPerNumber);
    for (T x : sample) {
        for (int b = 0; b < bitsPerNumber; ++b) {
            bits.push_back((x >> b) & 1);
        }
    }
    return bits;
}

/**
 * @brief Сохраняет сводку хи-квадрат тестов в CSV файл.
 * @tparam T Тип элементов выборки
 * @param filename Имя файла
 * @param name Имя генератора
 * @param samples Вектор выборок
 * @param minVal Минимальное значение генератора
 * @param maxVal Максимальное значение генератора
 * @param significance Уровень значимости
 */
template <typename T>
void saveChiSquareSummary(const std::string& filename,
                          const std::string& name,
                          const std::vector<std::vector<T>>& samples,
                          T minVal, T maxVal, double significance) {
    std::ofstream file(filename, std::ios::app); // append mode
    if (!file.is_open()) {
        // Если файл не открылся, пробуем создать с заголовком
        file.open(filename);
        file << "generator,passed_samples,total_samples,pass_rate\n";
    }
    
    size_t passed = 0;
    for (const auto& sample : samples) {
        bool chi_ok = custom::chi_squared_test(sample, minVal, maxVal,
                                               significance, custom::DistributionType::UNIFORM_INT);
        if (chi_ok) ++passed;
    }
    double passRate = static_cast<double>(passed) / samples.size();
    file << name << "," << passed << "," << samples.size() << "," << passRate << "\n";
}

/**
 * @brief Запускает пять тестов NIST на объединённой битовой последовательности.
 * @tparam T Тип элементов выборки
 * @param name Имя генератора
 * @param samples Вектор выборок
 * @param significance Уровень значимости
 */
template <typename T>
void runNISTTests(const std::string& name,
                  const std::vector<std::vector<T>>& samples,
                  double significance) {
    std::cout << "\nNIST tests for " << name << ":\n";

    std::vector<bool> allBits;
    for (const auto& sample : samples) {
        auto bits = toBits(sample, sizeof(T) * 8);
        allBits.insert(allBits.end(), bits.begin(), bits.end());
    }

    nist::FrequencyTest freq(allBits);
    nist::BlockFrequencyTest blockFreq(allBits, 8);
    nist::RunsTest runs(allBits);
    nist::LongestRunTest longestRun(allBits);
    nist::CumulativeSumsTest cumSums(allBits);

    auto status = [&](double p) { return (p >= significance) ? "PASS" : "FAIL"; };

    std::cout << "  Frequency test:        " << status(freq.p_value) << "\n";
    std::cout << "  Block frequency test:  " << status(blockFreq.p_value) << "\n";
    std::cout << "  Runs test:             " << status(runs.p_value) << "\n";
    std::cout << "  Longest run test:      " << status(longestRun.p_value) << "\n";
    std::cout << "  Cumulative sums test:  " << status(cumSums.p_value) << "\n";
}

/**
 * @brief Замеряет время генерации заданного количества чисел.
 * @tparam Generator Тип генератора
 * @param gen Генератор (передаётся по ссылке)
 * @param sizes Вектор объёмов выборок
 * @return Вектор времён в наносекундах
 */
template <typename Generator>
std::vector<uint64_t> measureTimes(Generator& gen, const std::vector<size_t>& sizes) {
    std::vector<uint64_t> times;
    for (size_t n : sizes) {
        auto start = std::chrono::high_resolution_clock::now();
        volatile typename Generator::result_type tmp;
        for (size_t j = 0; j < n; ++j) {
            tmp = gen.next();
        }
        auto end = std::chrono::high_resolution_clock::now();
        times.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
    }
    return times;
}

/**
 * @brief Сохраняет результаты замеров времени для модифицированных + mt19937_64.
 * @param filename Имя файла
 * @param sizes Объёмы выборок
 * @param t1 MLCG_XOR
 * @param t2 LFG_XOR
 * @param t3 Xorshift128_MLT
 * @param t4 mt19937_64
 */
void saveTimingsModified(const std::string& filename,
                         const std::vector<size_t>& sizes,
                         const std::vector<uint64_t>& t1,
                         const std::vector<uint64_t>& t2,
                         const std::vector<uint64_t>& t3,
                         const std::vector<uint64_t>& t4) {
    std::ofstream file(filename);
    if (!file.is_open()) return;
    file << "size,MLCG_XOR,LFG_XOR,Xorshift128_MLT,mt19937_64\n";
    for (size_t i = 0; i < sizes.size(); ++i) {
        file << sizes[i] << "," << t1[i] << "," << t2[i] << "," << t3[i] << "," << t4[i] << "\n";
    }
}

/**
 * @brief Сохраняет результаты замеров времени для всех генераторов (модифицированные + стандартные).
 * @param filename Имя файла
 * @param sizes Объёмы выборок
 * @param timings1..6 Векторы времён
 */
void saveTimingsComparison(const std::string& filename,
                           const std::vector<size_t>& sizes,
                           const std::vector<uint64_t>& t1,
                           const std::vector<uint64_t>& t2,
                           const std::vector<uint64_t>& t3,
                           const std::vector<uint64_t>& t4,
                           const std::vector<uint64_t>& t5,
                           const std::vector<uint64_t>& t6) {
    std::ofstream file(filename);
    if (!file.is_open()) return;
    file << "size,MLCG_XOR,LFG_XOR,Xorshift128_MLT,LCG_Standart,LFG_Standart,Xorshift128_Standart\n";
    for (size_t i = 0; i < sizes.size(); ++i) {
        file << sizes[i] << ","
             << t1[i] << "," << t2[i] << "," << t3[i] << ","
             << t4[i] << "," << t5[i] << "," << t6[i] << "\n";
    }
}

/**
 * @brief Сохраняет пары последовательных чисел для построения графика корреляций.
 * @tparam Generator Тип генератора
 * @param filename Имя файла
 * @param gen Генератор
 * @param numPoints Количество точек
 */
template <typename Generator>
void saveCorrelationData(const std::string& filename, Generator& gen, size_t numPoints = 1000) {
    std::ofstream file(filename);
    if (!file.is_open()) return;
 
    double scale = 1.0 / (static_cast<double>(gen.max()) + 1.0);
    
    double prev = static_cast<double>(gen.next()) * scale;
    for (size_t i = 0; i < numPoints; ++i) {
        double curr = static_cast<double>(gen.next()) * scale;
        file << prev << "," << curr << "\n";
        prev = curr;
    }
}

/**
 * @brief Тестирует генератор с разными seed (включая плохие) и сохраняет результат хи-квадрат.
 * @tparam Generator Тип генератора
 * @param seeds Вектор seed-значений для тестирования
 * @param genName Имя генератора (для вывода и имени файла)
 * @param numSamples Количество выборок на каждый seed
 * @param sampleSize Размер каждой выборки
 * @param minVal Минимальное значение генератора
 * @param maxVal Максимальное значение генератора
 * @param significance Уровень значимости
 * @param filename CSV-файл для сохранения результатов
 */
template <typename Generator>
void testSeedDependency(const std::vector<uint64_t>& seeds,
                        const std::string& genName,
                        size_t numSamples,
                        size_t sampleSize,
                        typename Generator::result_type  minVal,
                        typename Generator::result_type  maxVal,
                        double significance,
                        const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) return;
    file << "seed,passed_samples,total_samples,pass_rate\n";
    
    std::cout << "\nSeed dependency test for " << genName << ":\n";
    
    for (uint64_t seed : seeds) {
        Generator gen(seed);
        auto samples = generateSamples(gen, numSamples, sampleSize);
        
        size_t passed = 0;
        for (const auto& sample : samples) {
            bool chi_ok = custom::chi_squared_test(sample, minVal, maxVal,
                                                   significance, custom::DistributionType::UNIFORM_INT);
            if (chi_ok) ++passed;
        }
        double passRate = static_cast<double>(passed) / numSamples;
        file << seed << "," << passed << "," << numSamples << "," << passRate << "\n";
        std::cout << "  seed=" << seed << ": passed " << passed << "/" << numSamples
                  << " (" << std::fixed << std::setprecision(2) << passRate * 100 << "%)\n";
    }
}

int main() {
    MLCG_XOR g1(123456789);
    LFG_XOR g2(123456789);
    Xorshift128_MLT g3(123456789);

    auto mod_samples1 = generateSamples(g1, NUM_SAMPLES, SAMPLE_SIZE);
    auto mod_samples2 = generateSamples(g2, NUM_SAMPLES, SAMPLE_SIZE);
    auto mod_samples3 = generateSamples(g3, NUM_SAMPLES, SAMPLE_SIZE);

    printStats("MLCG_XOR", mod_samples1, MLCG_XOR::min(), MLCG_XOR::max(), SIGNIFICANCE);
    printStats("LFG_XOR", mod_samples2, LFG_XOR::min(), LFG_XOR::max(), SIGNIFICANCE);
    printStats("Xorshift128_MLT", mod_samples3, Xorshift128_MLT::min(), Xorshift128_MLT::max(), SIGNIFICANCE);

    runNISTTests("MLCG_XOR", mod_samples1, SIGNIFICANCE);
    runNISTTests("LFG_XOR", mod_samples2, SIGNIFICANCE);
    runNISTTests("Xorshift128_MLT", mod_samples3, SIGNIFICANCE);
    
    LCG_Standart s1(123456789);
    LFG_Standart s2(123456789);
    Xorshift128_Standart s3(123456789);

    auto std_samples1 = generateSamples(s1, NUM_SAMPLES, SAMPLE_SIZE);
    auto std_samples2 = generateSamples(s2, NUM_SAMPLES, SAMPLE_SIZE);
    auto std_samples3 = generateSamples(s3, NUM_SAMPLES, SAMPLE_SIZE);

    printStats("LCG_Standart", std_samples1, LCG_Standart::min(), LCG_Standart::max(), SIGNIFICANCE);
    printStats("LFG_Standart", std_samples2, LFG_Standart::min(), LFG_Standart::max(), SIGNIFICANCE);
    printStats("Xorshift128_Standart", std_samples3, Xorshift128_Standart::min(), Xorshift128_Standart::max(), SIGNIFICANCE);

    runNISTTests("LCG_Standart", std_samples1, SIGNIFICANCE);
    runNISTTests("LFG_Standart", std_samples2, SIGNIFICANCE);
    runNISTTests("Xorshift128_Standart", std_samples3, SIGNIFICANCE);

    saveChiSquareSummary("chi_square_summary.csv", "MLCG_XOR", mod_samples1, 
                     MLCG_XOR::min(), MLCG_XOR::max(), SIGNIFICANCE);
    saveChiSquareSummary("chi_square_summary.csv", "LFG_XOR", mod_samples2, 
                        LFG_XOR::min(), LFG_XOR::max(), SIGNIFICANCE);
    saveChiSquareSummary("chi_square_summary.csv", "Xorshift128_MLT", mod_samples3, 
                        Xorshift128_MLT::min(), Xorshift128_MLT::max(), SIGNIFICANCE);
    saveChiSquareSummary("chi_square_summary.csv", "LCG_Standart", std_samples1, 
                        LCG_Standart::min(), LCG_Standart::max(), SIGNIFICANCE);
    saveChiSquareSummary("chi_square_summary.csv", "LFG_Standart", std_samples2, 
                        LFG_Standart::min(), LFG_Standart::max(), SIGNIFICANCE);
    saveChiSquareSummary("chi_square_summary.csv", "Xorshift128_Standart", std_samples3, 
                        Xorshift128_Standart::min(), Xorshift128_Standart::max(), SIGNIFICANCE);

    MLCG_XOR corr1(123456789);
    LFG_XOR corr2(123456789);
    Xorshift128_MLT corr3(123456789);
    LCG_Standart corr4(123456789);
    LFG_Standart corr5(123456789);
    Xorshift128_Standart corr6(123456789);

    saveCorrelationData("corr_mlcg_xor.csv", corr1);
    saveCorrelationData("corr_lfg_xor.csv", corr2);
    saveCorrelationData("corr_xorshift128_mlt.csv", corr3);
    saveCorrelationData("corr_lcg_standart.csv", corr4);
    saveCorrelationData("corr_lfg_standart.csv", corr5);
    saveCorrelationData("corr_xorshift128_standart.csv", corr6);

    std::vector<size_t> sizes = {1000, 50000, 100000, 500000, 1000000};

    MLCG_XOR t1(123456789);
    LFG_XOR t2(123456789);
    Xorshift128_MLT t3(123456789);
    auto mod_t1 = measureTimes(t1, sizes);
    auto mod_t2 = measureTimes(t2, sizes);
    auto mod_t3 = measureTimes(t3, sizes);

    LCG_Standart t4(123456789);
    LFG_Standart t5(123456789);
    Xorshift128_Standart t6(123456789);
    auto std_t1 = measureTimes(t4, sizes);
    auto std_t2 = measureTimes(t5, sizes);
    auto std_t3 = measureTimes(t6, sizes);

    std::vector<uint64_t> mt_t;
    for (size_t n : sizes) {
        std::mt19937_64 mt_local(123456789);
        auto start = std::chrono::high_resolution_clock::now();
        volatile uint64_t tmp;
        for (size_t j = 0; j < n; ++j) {
            tmp = mt_local();
        }
        auto end = std::chrono::high_resolution_clock::now();
        mt_t.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
    }

    saveTimingsModified("timings_modified.csv", sizes, mod_t1, mod_t2, mod_t3, mt_t);
    
    saveTimingsComparison("timings_comparison.csv", sizes,
                          mod_t1, mod_t2, mod_t3,
                          std_t1, std_t2, std_t3);

    std::vector<uint64_t> bad_seeds = {0, 1, 12345, 0xFFFFFFFFFFFFFFFFULL};
    
    testSeedDependency<MLCG_XOR>(bad_seeds, "MLCG_XOR", NUM_SAMPLES, SAMPLE_SIZE,
                                 MLCG_XOR::min(), MLCG_XOR::max(), SIGNIFICANCE,
                                 "seed_mlcg_xor.csv");
    testSeedDependency<LFG_XOR>(bad_seeds, "LFG_XOR", NUM_SAMPLES, SAMPLE_SIZE,
                                LFG_XOR::min(), LFG_XOR::max(), SIGNIFICANCE,
                                "seed_lfg_xor.csv");
    testSeedDependency<Xorshift128_MLT>(bad_seeds, "Xorshift128_MLT", NUM_SAMPLES, SAMPLE_SIZE,
                                        Xorshift128_MLT::min(), Xorshift128_MLT::max(), SIGNIFICANCE,
                                        "seed_xorshift128_mlt.csv");
    
    testSeedDependency<LCG_Standart>(bad_seeds, "LCG_Standart", NUM_SAMPLES, SAMPLE_SIZE,
                                     LCG_Standart::min(), LCG_Standart::max(), SIGNIFICANCE,
                                     "seed_lcg_standart.csv");
    testSeedDependency<LFG_Standart>(bad_seeds, "LFG_Standart", NUM_SAMPLES, SAMPLE_SIZE,
                                     LFG_Standart::min(), LFG_Standart::max(), SIGNIFICANCE,
                                     "seed_lfg_standart.csv");
    testSeedDependency<Xorshift128_Standart>(bad_seeds, "Xorshift128_Standart", NUM_SAMPLES, SAMPLE_SIZE,
                                             Xorshift128_Standart::min(), Xorshift128_Standart::max(), SIGNIFICANCE,
                                             "seed_xorshift128_standart.csv");

    return 0;
}