# Renderer

Welcome to my 3rd Year Computer Graphics Coursework project.

[![Renderer Demo](https://img.youtube.com/vi/70PwL_9v6Z8/0.jpg)](https://youtube.com/shorts/70PwL_9v6Z8?feature=share)

> **The video above is a demo showcasing some of the features implemented in the renderer.**

## Overview

**RedNoise** is a custom renderer written in C++ without relying on any external rendering libraries. This project
demonstrates my understanding of computer graphics concepts and showcases features typically found in rendering
pipelines.

> **Note**:
> - The renderer is not yet optimized, but I plan to implement **Hierarchical Bounding Volumes (HBV)** in the future to
    improve performance and efficiency.
> - Please refer to the `CMakeLists.txt` file for detailed instructions on how to build and run the code.

## Features

- **OBJ File Loading and Parsing**  
  Efficiently loads and processes 3D models in the OBJ format.

- **Rasterizer**  
  A triangle rasterization pipeline for real-time rendering.

- **Ray Tracer**  
  A feature-rich ray tracing implementation, including:
    - **Soft Shadows**: Realistic light scattering for smooth shadow edges.
    - **Mirrors**: Reflective surfaces for enhanced scene realism.
    - **Texture Mapping**: Applying textures to 3D models.
    - **Phong Shading**: Simulating realistic lighting effects on surfaces.

- **Wireframe Mode**  
  Visualize 3D models in a clean, skeletal representation for better analysis and debugging.

---

Thank you for checking out my renderer!
