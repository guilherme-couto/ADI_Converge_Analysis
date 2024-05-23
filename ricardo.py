import numpy as np
import matplotlib.pyplot as plt

def load_from_file(filename):
    # Carregar os dados do arquivo
    data = np.loadtxt(filename)

    # Extrair as colunas do arquivo
    linhas = data[:, 0].astype(int)
    colunas = data[:, 1].astype(int)
    valores = data[:, 2]

    # Calcular o número de linhas e colunas
    num_linhas = np.max(linhas)
    num_colunas = np.max(colunas)

    # Criar uma matriz vazia para armazenar os valores
    matriz = np.zeros((num_linhas, num_colunas))

    # Preencher a matriz com os valores do arquivo
    for linha, coluna, valor in zip(linhas, colunas, valores):
        matriz[linha - 1, coluna - 1] = valor

    return matriz

# Nome do arquivo exata
filename = f'ricardo-exact.txt'

# Carregar os dados do arquivo
matriz_exata = load_from_file(filename)

# Plotar a matriz com imshow
plt.figure()
plt.imshow(matriz_exata, cmap='viridis', vmin=-1.0, vmax=1.0)
plt.colorbar(label='Value')
plt.xticks([])
plt.yticks([])
plt.title(f'Exata Ricardo dt=0.025 dx=0.05')
plt.savefig(f'./exata-ricardo-0.025-0.05.png')
plt.close()

# Nome do arquivo resultado
filename = f'ricardo-result.txt'

# Carregar os dados do arquivo
matriz_result = load_from_file(filename)

# Plotar a matriz com imshow
plt.figure()
plt.imshow(matriz_result, cmap='viridis', vmin=-1.0, vmax=1.0)
plt.colorbar(label='Value')
plt.xticks([])
plt.yticks([])
plt.title(f'Resultado Ricardo dt=0.025 dx=0.05')
plt.savefig(f'./resultado-ricardo-0.025-0.05.png')
plt.close()

# Calcular a diferença entre as duas matrizes
diferenca = np.abs(matriz_result - matriz_exata)
plt.figure()
plt.imshow(diferenca, cmap='viridis')
plt.colorbar(label='Value')
plt.xticks([])
plt.yticks([])
plt.title(f'Erros Ricardo dt=0.025 dx=0.05 (max_error={diferenca.max()})')
plt.savefig(f'./errors-ricardo-0.025-0.05.png')
plt.close()

# Diff between gui last and ricardo result
gui_last = np.genfromtxt(f'./simulation-files/double/AFHN/ADI/last-0.02500000-0.050000.txt', dtype=float)
diff_results = np.abs(gui_last - matriz_result)
plt.figure()
plt.imshow(diff_results, cmap='viridis')
plt.colorbar(label='Value')
plt.xticks([])
plt.yticks([])
plt.title(f'Diff Result Gui x Ricardo dt=0.025 dx=0.05 (max_error={diff_results.max()})')
plt.savefig(f'./diff-result-gui-ric-0.025-0.05.png')
plt.close()


# Read data from the text file
gui_exact = np.genfromtxt(f'./simulation-files/double/AFHN/ADI/exact-0.02500000-0.050000.txt', dtype=float)
diff_exacts = np.abs(gui_exact - matriz_exata)
plt.figure()
plt.imshow(diff_exacts, cmap='viridis')
plt.colorbar(label='Value')
plt.xticks([])
plt.yticks([])
plt.title(f'Diff Exacts Gui x Ricardo dt=0.025 dx=0.05 (max_error={diff_exacts.max()})')
plt.savefig(f'./diff-exact-gui-ric-0.025-0.05.png')
plt.close()

