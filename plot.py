import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Настройки стиля
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 10

# График 1: производительность (модифицированные + mt19937_64) 
def plot_performance_modified(csv_file='timings_modified.csv', output_file='performance_modified.png'):
    """
    Строит график производительности для модифицированных генераторов и std::mt19937_64.
    """
    df = pd.read_csv(csv_file)
    
    plt.figure(figsize=(10, 7))
    
    colors = {'MLCG_XOR': '#1f77b4', 'LFG_XOR': '#2ca02c', 
              'Xorshift128_MLT': '#d62728', 'mt19937_64': '#ff7f0e'}
    markers = {'MLCG_XOR': 'o', 'LFG_XOR': 's', 'Xorshift128_MLT': '^', 'mt19937_64': 'D'}
    
    for col in ['MLCG_XOR', 'LFG_XOR', 'Xorshift128_MLT', 'mt19937_64']:
        plt.plot(df['size'], df[col] / 1e6, 
                 marker=markers[col], label=col, linewidth=2, 
                 color=colors[col], markersize=6)
    
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Number of generated numbers (log scale)')
    plt.ylabel('Time (milliseconds, log scale)')
    plt.title('Performance comparison of modified generators vs std::mt19937_64')
    plt.legend(loc='upper left')
    plt.grid(True, which='both', linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    plt.show()
    print(f"Graph saved: {output_file}")

# График 2: сравнение модифицированных и стандартных 
def plot_performance_comparison(csv_file='timings_comparison.csv', output_file='performance_comparison.png'):
    """
    Строит график сравнения производительности всех 6 генераторов.
    Каждый генератор – свой цвет.
    """
    df = pd.read_csv(csv_file)
    
    plt.figure(figsize=(12, 7))
    
    # Цвета для каждого генератора (6 разных)
    colors = {
        'MLCG_XOR': '#1f77b4',      # синий
        'LFG_XOR': '#2ca02c',       # зелёный
        'Xorshift128_MLT': '#d62728', # красный
        'LCG_Standart': '#ff7f0e',   # оранжевый
        'LFG_Standart': '#9467bd',   # фиолетовый
        'Xorshift128_Standart': '#8c564b' # коричневый
    }
    
    markers = {
        'MLCG_XOR': 'o',
        'LFG_XOR': 's',
        'Xorshift128_MLT': '^',
        'LCG_Standart': 'D',
        'LFG_Standart': 'v',
        'Xorshift128_Standart': '<'
    }
    
    for col in df.columns:
        if col == 'size':
            continue
        plt.plot(df['size'], df[col] / 1e6, 
                 marker=markers[col], label=col, linewidth=2, 
                 color=colors[col], markersize=6)
    
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Number of generated numbers (log scale)')
    plt.ylabel('Time (milliseconds, log scale)')
    plt.title('Performance comparison of all generators')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9)
    plt.grid(True, which='both', linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    plt.show()
    print(f"Graph saved: {output_file}")

# График 3: корреляции с зумом 
def plot_correlations(csv_files, names, output_file='correlations.png'):
    """
    Строит 6 графиков корреляций (x_i, x_{i+1}) для всех генераторов.
    Для каждого генератора показывает полный диапазон и зум (0-0.1).
    """
    fig, axes = plt.subplots(2, 6, figsize=(18, 6))
    
    for idx, (csv_file, name) in enumerate(zip(csv_files, names)):
        df = pd.read_csv(csv_file, header=None, names=['x', 'y'])
        
        # Верхний ряд: полный диапазон [0, 1]
        ax_top = axes[0, idx]
        ax_top.scatter(df['x'], df['y'], s=0.5, alpha=0.3, c='blue')
        ax_top.set_title(f'{name}\n[0, 1]', fontsize=10)
        ax_top.set_xlabel('x_i')
        ax_top.set_ylabel('x_{i+1}')
        ax_top.set_xlim(0, 1)
        ax_top.set_ylim(0, 1)
        ax_top.set_aspect('equal')
        
        # Нижний ряд: зум [0, 0.1] для выявления структур
        ax_bottom = axes[1, idx]
        ax_bottom.scatter(df['x'], df['y'], s=1, alpha=0.5, c='blue')
        ax_bottom.set_title(f'{name}\nzoom [0, 0.1]', fontsize=10)
        ax_bottom.set_xlabel('x_i')
        ax_bottom.set_ylabel('x_{i+1}')
        ax_bottom.set_xlim(0, 0.1)
        ax_bottom.set_ylim(0, 0.1)
        ax_bottom.set_aspect('equal')
        ax_bottom.grid(True, linestyle='--', alpha=0.5)
    
    plt.suptitle('Correlation plots: (x_i, x_{i+1})', fontsize=14)
    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    plt.show()
    print(f"Graph saved: {output_file}")

# График 4: зависимость от seed 
def plot_seed_dependency(seed_files, names, output_file='seed_dependency.png'):
    """
    Строит гистограммы прохождения хи-квадрат теста для разных seed.
    """
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()
    
    for idx, (csv_file, name) in enumerate(zip(seed_files, names)):
        df = pd.read_csv(csv_file)
        ax = axes[idx]
        
        seed_labels = [str(s) for s in df['seed']]
        seed_labels = ['2^64-1' if s == '18446744073709551615' else s for s in seed_labels]
        
        bars = ax.bar(seed_labels, df['pass_rate'], color='steelblue', edgecolor='black')
        ax.set_ylim(0, 1.05)
        ax.set_ylabel('Pass rate (Chi-square)')
        ax.set_title(f'{name}')
        ax.axhline(y=0.95, color='red', linestyle='--', linewidth=2, label='expected (95%)')
        ax.legend()
        
        # Подписи значений над столбцами
        for bar, rate in zip(bars, df['pass_rate']):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                    f'{rate:.2f}', ha='center', va='bottom', fontsize=9)
        
        ax.tick_params(axis='x', rotation=15)
    
    for idx in range(len(seed_files), 6):
        axes[idx].set_visible(False)
    
    plt.suptitle('Seed dependency: Chi-square test pass rate', fontsize=14)
    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    plt.show()
    print(f"Graph saved: {output_file}")

# График 5: статистика прохождения хи-квадрат по выборкам 
def plot_chi_square_summary(csv_file='chi_square_summary.csv', output_file='chi_square_summary.png'):
    """
    Строит сводную диаграмму доли выборок, прошедших хи-квадрат тест.
    Данные читаются из CSV, созданного C++ программой.
    """
    try:
        df = pd.read_csv(csv_file)
    except FileNotFoundError:
        print(f"Warning: {csv_file} not found. Run C++ program first to generate data.")
        return
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Сортируем для красоты
    df = df.sort_values('pass_rate', ascending=False)
    
    colors = plt.cm.Set3(range(len(df)))
    bars = ax.bar(df['generator'], df['pass_rate'], color=colors, edgecolor='black')
    
    ax.set_ylim(0, 1.05)
    ax.set_ylabel('Proportion of samples passing Chi-square test')
    ax.set_title('Chi-square test pass rate by generator (α=0.05)')
    
    # Подписи значений над столбцами
    for bar, rate in zip(bars, df['pass_rate']):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                f'{rate:.2f}', ha='center', va='bottom', fontsize=10)
    
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    plt.show()
    print(f"Graph saved: {output_file}")

# График 6: сравнение корреляций LCG_Standart и MLCG_XOR 
def plot_correlation_comparison(csv_file_bad='corr_lcg_standart.csv', 
                                 csv_file_good='corr_mlcg_xor.csv',
                                 output_file='correlation_comparison.png'):
    """
    Строит два графика корреляций рядом: для LCG_Standart и MLCG_XOR.
    Добавляет область с зумом (0-0.1) и статистический анализ.
    """
    df_bad = pd.read_csv(csv_file_bad, header=None, names=['x', 'y'])
    df_good = pd.read_csv(csv_file_good, header=None, names=['x', 'y'])
    
    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)
    
    # LCG_Standart (полный диапазон)
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.scatter(df_bad['x'], df_bad['y'], s=0.5, alpha=0.3, c='red')
    ax1.set_title('LCG_Standart (full range [0,1])', fontsize=12)
    ax1.set_xlabel('x_i')
    ax1.set_ylabel('x_{i+1}')
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1)
    ax1.set_aspect('equal')
    
    # MLCG_XOR (полный диапазон)
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.scatter(df_good['x'], df_good['y'], s=0.5, alpha=0.3, c='blue')
    ax2.set_title('MLCG_XOR (full range [0,1])', fontsize=12)
    ax2.set_xlabel('x_i')
    ax2.set_ylabel('x_{i+1}')
    ax2.set_xlim(0, 1)
    ax2.set_ylim(0, 1)
    ax2.set_aspect('equal')
    
    # LCG_Standart (zoom)
    ax3 = fig.add_subplot(gs[1, 0])
    ax3.scatter(df_bad['x'], df_bad['y'], s=1, alpha=0.5, c='red')
    ax3.set_title('LCG_Standart (zoom [0,0.1])', fontsize=12)
    ax3.set_xlabel('x_i')
    ax3.set_ylabel('x_{i+1}')
    ax3.set_xlim(0, 0.1)
    ax3.set_ylim(0, 0.1)
    ax3.set_aspect('equal')
    ax3.grid(True, linestyle='--', alpha=0.5)
    
    # MLCG_XOR (zoom)
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.scatter(df_good['x'], df_good['y'], s=1, alpha=0.5, c='blue')
    ax4.set_title('MLCG_XOR (zoom [0,0.1])', fontsize=12)
    ax4.set_xlabel('x_i')
    ax4.set_ylabel('x_{i+1}')
    ax4.set_xlim(0, 0.1)
    ax4.set_ylim(0, 0.1)
    ax4.set_aspect('equal')
    ax4.grid(True, linestyle='--', alpha=0.5)
    
    # Статистический анализ 
    def calculate_stats(df, name):
        # Корреляция Пирсона
        corr = df['x'].corr(df['y'])
        # Коэффициент детерминации
        r2 = corr ** 2
        return f"{name}\nPearson r = {corr:.4f}\nR² = {r2:.4f}"
    
    stats_bad = calculate_stats(df_bad, "LCG_Standart")
    stats_good = calculate_stats(df_good, "MLCG_XOR")
    
    ax1.text(0.05, 0.95, stats_bad, transform=ax1.transAxes, fontsize=9,
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax2.text(0.05, 0.95, stats_good, transform=ax2.transAxes, fontsize=9,
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.suptitle('Correlation comparison: LCG_Standart (unmodified) vs MLCG_XOR (modified)', fontsize=14)
    plt.savefig(output_file, dpi=150)
    plt.show()
    print(f"Graph saved: {output_file}")

if __name__ == '__main__':
    print("Starting graph generation...")
    print("=" * 50)
    
    # 1. График производительности (модифицированные + mt19937_64)
    try:
        plot_performance_modified('timings_modified.csv', 'performance_modified.png')
    except FileNotFoundError:
        print("Warning: timings_modified.csv not found. Skipping performance_modified graph.")
    
    # 2. График сравнения модифицированных и стандартных
    try:
        plot_performance_comparison('timings_comparison.csv', 'performance_comparison.png')
    except FileNotFoundError:
        print("Warning: timings_comparison.csv not found. Skipping performance_comparison graph.")
    
    # 3. Графики корреляций (6 генераторов)
    corr_files = [
        'corr_mlcg_xor.csv',
        'corr_lfg_xor.csv',
        'corr_xorshift128_mlt.csv',
        'corr_lcg_standart.csv',
        'corr_lfg_standart.csv',
        'corr_xorshift128_standart.csv'
    ]
    corr_names = ['MLCG_XOR', 'LFG_XOR', 'Xorshift128_MLT', 
                  'LCG_Standart', 'LFG_Standart', 'Xorshift128_Standart']
    
    try:
        plot_correlations(corr_files, corr_names, 'correlations.png')
    except FileNotFoundError as e:
        print(f"Warning: {e}. Skipping correlations graph.")
    
    # 4. Графики зависимости от seed
    seed_files = [
        'seed_mlcg_xor.csv',
        'seed_lfg_xor.csv',
        'seed_xorshift128_mlt.csv',
        'seed_lcg_standart.csv',
        'seed_lfg_standart.csv',
        'seed_xorshift128_standart.csv'
    ]
    seed_names = ['MLCG_XOR', 'LFG_XOR', 'Xorshift128_MLT',
                  'LCG_Standart', 'LFG_Standart', 'Xorshift128_Standart']
    
    try:
        plot_seed_dependency(seed_files, seed_names, 'seed_dependency.png')
    except FileNotFoundError as e:
        print(f"Warning: {e}. Skipping seed_dependency graph.")
    
    # 5. Сводная диаграмма хи-квадрат (требует ручного ввода данных)
    # plot_chi_square_summary([], [], 'chi_square_summary.png')
    
    # 6. Сравнение корреляций LCG_Standart и MLCG_XOR
    try:
        plot_correlation_comparison('corr_lcg_standart.csv', 'corr_mlcg_xor.csv', 
                                    'correlation_comparison.png')
    except FileNotFoundError as e:
        print(f"Warning: {e}. Skipping correlation_comparison graph.")
    
    print("=" * 50)
    print("Graph generation completed!")