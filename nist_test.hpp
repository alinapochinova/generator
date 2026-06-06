/**
 * @file nist_test.hpp
 * @brief Реализация пяти статистических тестов из набора NIST SP 800-22.
 * @details Набор тестов предназначен для проверки случайности битовых последовательностей.
 *          Реализованы следующие тесты:
 *          - Frequency (Monobit) Test          (проверка глобальной доли единиц)
 *          - Block Frequency Test              (проверка равномерности внутри блоков)
 *          - Runs Test                         (проверка количества серий)
 *          - Longest Run Test                  (проверка максимальной длины серии)
 *          - Cumulative Sums Test              (проверка случайности накопленной суммы)
 * 
 *          Все тесты вычисляют p-value. Если p-value ge заданного уровня значимости, 
 *          последовательность считается случайной.
 */

# pragma once

#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

namespace nist {

    /**
     * @brief Дополнительная функция ошибок (erfc)
     * @param x аргумент.
     * @return 1 - erf(x)
     */
    inline double erfc(double x) {
        return 1 - std::erf(x);
    }

    /**
     * @brief Частотный (monobit) тест — проверяет, близка ли доля единиц во всей
     *        последовательности к 1/2.
     * @details Статистика: S_n = (количество единиц) - (количество нулей).
     *          Нормированная статистика: S = |S_n| / sqrtn.
     *          p-value = erfc(S/sqrt2). Если p-value ge alpha, последовательность случайна.
     */
    struct FrequencyTest {
        double p_value;

        /**
         * @param bits входная битовая последовательность.
         */
        FrequencyTest(const std::vector<bool>& bits) {
            int64_t sum = 0;
            for (bool b : bits) sum += (b ? 1 : -1);
            double S = std::abs(sum) / std::sqrt(bits.size());
            p_value = erfc(S / std::sqrt(2.0));
        }
    };

    /**
     * @brief Блочный частотный тест — проверяет равномерность распределения единиц
     *        внутри блоков фиксированного размера.
     * @details Последовательность разбивается на N непересекающихся блоков размера M.
     *          Для каждого блока вычисляется доля единиц pi_i. Статистика:
     *          chi-squared = 4M sum (pi_i - 0.5)². Приближённо распределена как chi-squared с N степенями
     *          свободы. Для больших N используется нормальная аппроксимация.
     */

    struct BlockFrequencyTest {
        double p_value;

        /**
         * @param bits входная последовательность.
         * @param blockSize размер блока (в битах).
         */
        BlockFrequencyTest(const std::vector<bool>& bits, size_t blockSize) {
            size_t M = blockSize;
            size_t N = bits.size() / M;
            if (N == 0) { p_value = 0; return; }
            double chi2 = 0.0;
            // Цикл проходит по каждому блоку по порядку
            for (size_t i = 0; i < N; ++i) {
                size_t cnt = 0; // счётчик единиц в текущем блоке
                // Внутренний цикл перебирает M битов текущего блока
                for (size_t j = 0; j < M; ++j)
                // Если бит равен 1 (истина), то cnt увеличивается
                    if (bits[i * M + j]) ++cnt;
                double pi = static_cast<double>(cnt) / M; // выборочная доля единиц в текущем блоке
                chi2 += (pi - 0.5) * (pi - 0.5);
            }
            chi2 *= 4.0 * M;
            // Для больших N используем нормальную аппроксимацию
            double mean = N;
            double var = 2 * N;
            double x = (chi2 - mean) / std::sqrt(var);
            p_value = erfc(std::abs(x) / std::sqrt(2.0));
        }
    };

    /**
     * @brief Тест на серии (runs test) — проверяет, соответствует ли количество серий
     *        (последовательностей одинаковых битов) случайному распределению.
     * @details Пусть pi = количество единиц / n. Если pi не слишком отличается от 0.5,
     *          статистика: Z = (R - 2npi(1-pi)) / (2sqrt(2n)pi(1-pi)), где R — общее число серий.
     *          p-value = erfc(|Z|/sqrt2). Малое p-value указывает на корреляцию соседних битов.
     */
    
    struct RunsTest {
        double p_value;

        RunsTest(const std::vector<bool>& bits) {
            // Подсчитываем количество единиц во всей последовательности
            int64_t ones = std::count(bits.begin(), bits.end(), true);
            // Вычисляем выборочную долю единиц pi
            double pi = static_cast<double>(ones) / bits.size();

            // Предварительная проверка: если pi сильно отклоняется от 0.5, тест не проводится
            // 2/sqrtn – это эмпирическое правило, используемое в NIST
            if (std::abs(pi - 0.5) > 2.0 / std::sqrt(bits.size())) { 
                p_value = 0.0;
                return;
            }

            // Начинаем подсчёт серий. Первая серия всегда существует, поэтому runs начинается с 1
            int64_t runs = 1;
            // Проходим по последовательности, начиная со второго бита
            for (size_t i = 1; i < bits.size(); ++i)
            // Если текущий бит отличается от предыдущего, значит, началась новая серия – увеличиваем счётчик runs
                if (bits[i] != bits[i-1]) ++runs;
            // Вычисляем числитель статистики
            double num = static_cast<double>(runs) - 2.0 * bits.size() * pi * (1.0 - pi);
            //Знаменатель
            double den = 2.0 * std::sqrt(2.0 * bits.size()) * pi * (1.0 - pi);
            p_value = erfc(std::abs(num) / den);
        }
    };

    /**
     * @brief Тест на самую длинную серию — проверяет, не является ли максимальная длина
     *        серии в блоках слишком большой (что указывает на неслучайность).
     * @details Последовательность разбивается на блоки по 8 бит. Для каждого блока
     *          находится максимальная длина серии (подряд идущих одинаковых битов).
     *          Затем сравнивается наблюдаемое распределение максимальных длин с теоретическим.
     *          Статистика chi-squared рассчитывается по 7 категориям (длина le1, 2, 3, 4, 5, 6, ge7).
     *          p-value вычисляется приближённо как exp(-chi-squared/2).
     */
    struct LongestRunTest {
        double p_value;

        LongestRunTest(const std::vector<bool>& bits) {
            size_t blockSize = 8;
            size_t N = bits.size() / blockSize;
            std::vector<int> maxRuns(N, 0); // вектор для хранения максимальной длины серии единиц (run) в каждом блоке

            // Для каждого блока i перебираются биты
            for (size_t i = 0; i < N; ++i) {
                int cur = 0, best = 0;
                for (size_t j = 0; j < blockSize; ++j) {
                    // Если бит равен 1, текущая длина cur увеличивается; если 0 – обнуляется
                    if (bits[i*blockSize + j]) cur++;
                    else cur = 0;
                    //В переменной best сохраняется максимальная длина серии единиц в этом блоке
                    if (cur > best) best = cur;
                }
                maxRuns[i] = best;
            }
            // Категории: 0–1, 2, 3, 4, 5, 6, ge7
            std::vector<int> obs(7,0);
            for (int r : maxRuns) {
                if (r <= 1) obs[0]++;
                else if (r == 2) obs[1]++;
                else if (r == 3) obs[2]++;
                else if (r == 4) obs[3]++;
                else if (r == 5) obs[4]++;
                else if (r == 6) obs[5]++;
                else obs[6]++;
            }
            // Теоретические вероятности для блоков длиной 8 бит (взяты из документа NIST SP 800‑22)
            std::vector<double> exp = {0.2148, 0.3672, 0.2305, 0.0977, 0.0312, 0.0084, 0.0502};
            double chi2 = 0.0;
            for (size_t i = 0; i < 7; ++i) {
                double e = exp[i] * N; // ожидаемое количество блоков в i‑й категории
                // obs[i] – наблюдаемое
                chi2 += (obs[i] - e) * (obs[i] - e) / e;
            }
            // Приближённое p-value
            p_value = std::exp(-chi2 / 2);
        }
    };

    /**
     * @brief Тест кумулятивных сумм — проверяет, насколько сильно отклоняется
     *        накопленная сумма (присваивая биту '1' +1, биту '0' –1) от нуля.
     * @details Вычисляются максимальное положительное и отрицательное отклонения
     *          S_max и S_min. Нормированные значения: z1 = S_max / sqrtn, z2 = |S_min| / sqrtn.
     *          p-value = P(z) = sum_{k=-infinity}^{infinity} [erfc((2k+1)z/sqrt2) - erfc((2k-1)z/sqrt2)].
     *          Окончательный p-value = min(P(z1), P(z2)). Малое p-value означает,
     *          что накопленная сумма слишком часто (или слишком редко) пересекает нуль.
     */
    struct CumulativeSumsTest {
        double p_value;

        CumulativeSumsTest(const std::vector<bool>& bits) {
            int64_t S = 0; // текущая кумулятивная сумма
            int64_t sup = 0, inf = 0;
            // Цикл проходит по всем битам последовательности
            for (bool b : bits) {
                // если бит = 1, прибавляем +1; если 0, вычитаем 1
                S += (b ? 1 : -1);
                // Обновляем sup и inf, если текущее значение вышло за предыдущие границы
                if (S > sup) sup = S;
                if (S < inf) inf = S;
            }
            // Нормируем отклонения на квадратный корень из длины последовательности
            double z1 = std::abs(sup) / std::sqrt(bits.size()); // нормированный положительный максимум
            double z2 = std::abs(inf) / std::sqrt(bits.size()); // нормированный отрицательный минимум

            // P(z) – это сумма по нечётным членам, дающая вероятность того, что максимум (по модулю) меньше z
            auto P = [](double z) {
                double sum = 0;
                // Суммирование по k от -5 до 5 даёт достаточную точность, так как члены с большими |k| пренебрежимо малы
                for (int k = -5; k <= 5; ++k) {
                    // разность дополнительных функций ошибок, которая эквивалентна разности нормальных CDF
                    double term = erfc((2*k+1)*z / std::sqrt(2.0)) - erfc((2*k-1)*z / std::sqrt(2.0));
                    sum += term;
                }
                return sum;
            };
            // итоговый p‑value – это вероятность того, что при случайном блуждании и положительное, и отрицательное отклонения будут не больше наблюдаемых
            p_value = std::min(P(z1), P(z2));
        }
    };

}
