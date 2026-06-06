// standart_generators.hpp
#ifndef STANDART_GENERATORS_HPP
#define STANDART_GENERATORS_HPP

#include <cstdint>

/**
 * @file standart_generators.hpp
 * @brief Немодифицированные (классические) версии генераторов для сравнения.
 * 
 * Содержит классические реализации LCG, LFG и Xorshift128
 * без каких-либо нелинейных модификаций (перемешивание, умножение и т.д.).
 * 
 * Цель: сравнить статистические свойства и производительность
 * модифицированных и немодифицированных версий.
 */

/**
 * @brief Классический линейный конгруэнтный генератор (LCG) без модификаций.
 * 
 * Формула: state = A * state + C mod 2^64.
 * Параметры A и C выбраны так же, как в модифицированной версии,
 * но без пост-перемешивания (XOR со сдвигами).
 * 
 * @note Период: 2^64 при правильных параметрах.
 * @warning Младшие биты имеют короткий период, есть корреляции между соседними числами.
 * 
 * @see MLCG_XOR – модифицированная версия с перемешиванием
 */
class LCG_Standart {
public:
    using result_type = uint64_t;

    /**
     * @brief Конструктор.
     * @param seed Начальное значение (по умолчанию 123456789)
     */
    explicit LCG_Standart(uint64_t seed = 123456789);

    /** @brief Генерирует следующее 64-битное псевдослучайное число. */
    uint64_t next();

    static constexpr uint64_t min() { return 0; }
    static constexpr uint64_t max() { return ~0ULL; }

private:
    uint64_t state;
    static constexpr uint64_t A = 0x9e3779b97f4a7c15ULL;
    static constexpr uint64_t C = 0xbf58476d1ce4e5b9ULL;
};

/**
 * @brief Классический генератор Фибоначчи с запаздыванием (LFG) без модификаций.
 * 
 * Формула: buffer[i] = buffer[i] ^ buffer[(i+24)%55] (без scramble).
 * Лаги: p=55, q=24.
 * 
 * @note Период: ~2^79 для XOR-варианта.
 * @warning Сильная зависимость от начального seed, возможны корреляции.
 * 
 * @see LFG_XOR – модифицированная версия с scramble и улучшенной инициализацией
 */
class LFG_Standart {
public:
    using result_type = uint32_t;

    /**
     * @brief Конструктор.
     * @param seed Начальное значение (по умолчанию 123456789)
     */
    explicit LFG_Standart(uint32_t seed = 123456789);

    /** @brief Генерирует следующее 32-битное псевдослучайное число. */
    uint32_t next();

    static constexpr uint32_t min() { return 0; }
    static constexpr uint32_t max() { return ~0U; }

private:
    static constexpr int LAG1 = 55;
    static constexpr int LAG2 = 24;
    uint32_t buffer[LAG1];
    int index;

    /** @brief Инициализация буфера с заданным seed (простая версия без scramble). */
    void initBuffer(uint32_t seed);
};

/**
 * @brief Классический Xorshift128 без нелинейного умножения.
 * 
 * Формула: только XOR и сдвиги (23, 17, 26).
 * Состояние: 128 бит (два 64-битных слова).
 * 
 * @note Период: 2^128 - 1.
 * @warning Линейность в GF(2) – по нескольким выходам можно восстановить состояние.
 * 
 * @see Xorshift128_MLT – модифицированная версия с умножением (SplitMix64)
 */
class Xorshift128_Standart {
public:
    using result_type = uint64_t;

    /**
     * @brief Конструктор.
     * @param seed Начальное значение (по умолчанию 123456789)
     */
    explicit Xorshift128_Standart(uint64_t seed = 123456789);

    /** @brief Генерирует следующее 64-битное псевдослучайное число. */
    uint64_t next();

    static constexpr uint64_t min() { return 0; }
    static constexpr uint64_t max() { return ~0ULL; }

private:
    uint64_t s[2];

    /** @brief Прогрев генератора (10 итераций для честного сравнения с модифицированной версией). */
    void warmup(int iterations = 10);
};

#endif // STANDART_GENERATORS_HPP