---
title: 'Curso de Genômica Humana'
date: 2025
type: landing

design:
  # Default section spacing
  spacing: "6rem"

sections:
  - block: hero
    content:
      title: Curso de Genômica Humana
      text: Fundamentos da análise de dados de sequenciamento de exoma completo
      primary_action:
        text: Comece aqui
        url: /CursoWesICC/docs/
        icon: 🧬
      secondary_action:
        text: Sobre os autores
        url: /CursoWesICC/showcase/
      announcement:
        text: "Para download dos datasets de estudo."
        link:
          text: "Clique aqui"
          url: "/CursoWesICC/blog/"
    design:
      spacing:
        padding: [0, 0, 0, 0]
        margin: [0, 0, 0, 0]
      # For full-screen, add `min-h-screen` below
      css_class: ""
      background:
        color: ""
        image:
          # Add your image background to `assets/media/`.
          filename: ""
          filters:
            brightness: 0.5
  - block: features
    id: features
    content:
      title: Conteúdo
      text: Conteúdo programático abordado no curso. (v2025-07.1)
      items:
        - name: Estrutura FASTQ
          icon: magnifying-glass
          description: Estrutura do arquivo de leituras (FASTQ).
        - name: QC pré-mapeamento
          icon: bolt
          description: Controle de qualidade de arquivos de leituras (FASTQ) nos programas FastQC/MultiQC.
        - name: Estrutura VCF
          icon: sparkles
          description: Estrutura do arquivo de chamada de variantes genéticas (VCF).
        - name: Exploração de dados em VCF 
          icon: code-bracket
          description: Interpretação das informações contidas em um VCF anotado.
        - name: Visualização de dados em VCF
          icon: star
          description: Aplicação de pacotes de visualização de dados no R para análise de variantes genéticas.
        - name: Exploração genômica orientada por hipóteses
          icon: rectangle-group
          description: Delineamento experimental com base em hipóteses biológicas.
---
