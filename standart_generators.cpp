// standart_generators.cpp
#include "standart_generators.hpp"

/**
 * @brief Конструктор LCG_Standart.
 * @param seed Начальное значение состояния.
 * 
 * Инициализирует внутреннее состояние заданным seed.
 * Состояние будет обновляться по формуле LCG при каждом вызове next().
 */
LCG_Standart::LCG_Standart(uint64_t seed) : state(seed) {}

/**
 * @brief Генерация следующего 64-битного числа.
 * @return Очередное значение LCG без пост-обработки.
 * 
 * Выполняется классический LCG шаг: state = A * state + C.
 * Переполнение uint64_t автоматически даёт модуль 2^64.
 * 
 * @note Отсутствие перемешивания приводит к корреляции младших битов.
 */
uint64_t LCG_Standart::next() {
    state = A * state + C;
    return state;
}

/**
 * @brief Конструктор LFG_Standart.
 * @param seed Начальное значение для инициализации буфера.
 * 
 * Устанавливает index = LAG1 (55), чтобы при первом вызове next()
 * буфер был перегенерирован. Вызывает initBuffer(seed) для заполнения
 * начальных значений буфера.
 */
LFG_Standart::LFG_Standart(uint32_t seed) : index(LAG1) {
    initBuffer(seed);
}

/**
 * @brief Простая инициализация буфера без перемешивания.
 * @param seed Начальное значение.
 * 
 * Заполняет буфер из 55 элементов с помощью LCG:
 * state = state * 1103515245 + 12345.
 * 
 * @note Нет функции scramble, нет дополнительных проходов перемешивания.
 * Это делает генератор чувствительным к выбору seed.
 */
void LFG_Standart::initBuffer(uint32_t seed) {
    uint32_t state = seed;
    for (int i = 0; i < LAG1; ++i) {
        state = state * 1103515245 + 12345;
        buffer[i] = state;   // Без scramble!
    }
    // Нет дополнительных проходов
}

/**
 * @brief Генерация следующего 32-битного числа.
 * @return Очередное значение LFG без scramble.
 * 
 * Если индекс достиг конца буфера (55), весь буфер пересчитывается
 * по формуле: buffer[i] ^= buffer[(i+24)%55].
 * Затем возвращается текущий элемент буфера, и индекс увеличивается.
 * 
 * @note Отсутствие scramble может приводить к корреляциям.
 */
uint32_t LFG_Standart::next() {
    if (index >= LAG1) {
        for (int i = 0; i < LAG1; ++i) {
            int j = (i + LAG2) % LAG1;
            buffer[i] ^= buffer[j];   // Только XOR, без scramble!
        }
        index = 0;
    }
    return buffer[index++];
}


/**
 * @brief Конструктор Xorshift128_Standart.
 * @param seed Начальное значение.
 * 
 * Инициализирует состояние: s[0] = seed, s[1] = seed * золотое сечение.
 * Затем выполняет прогрев (10 итераций) для честного сравнения
 * с модифицированной версией Xorshift128_MLT.
 */
Xorshift128_Standart::Xorshift128_Standart(uint64_t seed) {
    s[0] = seed;
    s[1] = seed * 0x9e3779b97f4a7c15ULL;
    warmup(10);
}

/**
 * @brief Прогрев генератора.
 * @param iterations Количество итераций прогрева.
 * 
 * Выполняет несколько итераций next() без использования результатов.
 * Это улучшает начальное состояние, особенно при плохих seed.
 */
void Xorshift128_Standart::warmup(int iterations) {
    for (int i = 0; i < iterations; ++i) {
        next();
    }
}

/**
 * @brief Генерация следующего 64-битного числа.
 * @return Очередное значение Xorshift без нелинейного умножения.
 * 
 * Выполняется классический Xorshift128 с параметрами (23, 17, 26):
 * - s1 ^= s1 << 23
 * - s1 ^= s1 >> 17
 * - s1 ^= s0
 * - s1 ^= s0 >> 26
 * 
 * @note Операции линейны в GF(2), что может приводить к корреляциям.
 * Модифицированная версия добавляет умножение для нелинейности.
 */
uint64_t Xorshift128_Standart::next() {
    uint64_t s1 = s[0];
    uint64_t s0 = s[1];
    s[0] = s0;
    s1 ^= s1 << 23;
    s1 ^= s1 >> 17;
    s1 ^= s0;
    s1 ^= s0 >> 26;
    s[1] = s1;
    return s1;   // Без умножения
}