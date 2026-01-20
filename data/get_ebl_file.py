import os
import shutil
import sys
import ebltable.ebl_from_model as ebl

# 1. Localizar onde o ebltable foi instalado
package_dir = os.path.dirname(ebl.__file__)
data_dir_source = os.path.join(package_dir, 'data')

# Nome do arquivo alvo
filename = 'tau_dominguez11.fits'
source_file = os.path.join(data_dir_source, filename)

# Destino no seu projeto
dest_dir = "data/ebl_models"
dest_file = os.path.join(dest_dir, filename)

print(f"--- Extrator de Arquivos EBL ---")
print(f"Procurando em: {source_file}")

if os.path.exists(source_file):
    # Criar pasta de destino se não existir
    if not os.path.exists(dest_dir):
        os.makedirs(dest_dir)
        
    # Copiar
    shutil.copy2(source_file, dest_file)
    print(f"\n[SUCESSO] Arquivo copiado para: {dest_file}")
    print("Agora você pode rodar o utils/ebl.py normalmente!")
else:
    print(f"\n[ERRO] Arquivo não encontrado no pacote instalado.")
    print(f"Conteúdo da pasta data do pacote: {os.listdir(data_dir_source)}")
    
    # Plano B: Tentar encontrar qualquer arquivo .fits que pareça ser o do Dominguez
    print("\nTentando busca inteligente...")
    for f in os.listdir(data_dir_source):
        if 'dominguez' in f.lower() and 'tau' in f.lower():
            print(f"Achado alternativo: {f}")
            shutil.copy2(os.path.join(data_dir_source, f), os.path.join(dest_dir, f))
            print(f"Copiado como substituto. Atualize o utils/ebl.py com o nome: {f}")