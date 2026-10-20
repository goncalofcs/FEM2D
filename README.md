# INSTRUÇÕES

## Código python:

Para o correto funcionamento do programa são necessários 3 ficheiros .py:

1. reader.py - Classe que permite a leitura do ficheiro de texto com os dados do problema.
2. visualize.py - Classe que permite a visualização da estrutura deformada e original, semelhante à visualização no Abaqus.
3. FEM2D_main.py - Ficheiro que contém o programa em si que deve ser executado.

OBS:
A linha 134, pode sofrer alterações dependendo do caminho que o ficheiro de texto (data.txt) foi guardado. Em caso do caminho ser diferente, editar essa linha direcionando para o caminho correto: `reader = Reader(filename="data.txt")`

## Ficheiros:

S_Parte1.rpt e U_RF_Parte1.rpt - Correspondem aos valores obtidos pelo Abaqus para a primeira parte do projeto em unidades S.I.
U_RF_Parte2.rpt - Corresponde aos valores obtidos pelo Abaqus para a segunda pate do projeto (deslocamento imposto de -2 mm em D) em unidades S.I.
Projeto1 - Corresponde ao ficheiro .cae do Abaqus.
