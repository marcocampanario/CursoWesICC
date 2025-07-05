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
  - block: stats
    content:
      items:
        - statistic: "1M+"
          description: |
            Websites built  
            with Hugo Blox
        - statistic: "10k+"
          description: |
            GitHub stars  
            since 2016
        - statistic: "3k+"
          description: |
            Discord community  
            for support
    design:
      # Section background color (CSS class)
      css_class: "bg-gray-100 dark:bg-gray-800"
      # Reduce spacing
      spacing:
        padding: ["1rem", 0, "1rem", 0]
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
  - block: cta-card
    content:
      title: "Start Writing with the #1 Effortless Documentation Platform"
      text: Hugo Blox Docs Theme brings all your technical knowledge together in a single, centralized knowledge base. Easily search and edit it with the tools you use every day!
      button:
        text: Get Started
        url: https://hugoblox.com/templates/details/docs/
    design:
      card:
        # Card background color (CSS class)
        css_class: "bg-primary-700"
        css_style: ""
---
