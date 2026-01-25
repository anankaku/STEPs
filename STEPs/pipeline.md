```mermaid
flowchart TD
    classDef file fill:#f9f,stroke:#333,stroke-width:2px,color:black;
    classDef tool fill:#bbf,stroke:#333,stroke-width:2px,color:black;
    classDef process fill:#e1f5fe,stroke:#333,stroke-width:2px,color:black,stroke-dasharray: 5 5;
    classDef note fill:#fff3cd,stroke:#333,stroke-width:1px,color:black;

    %% =========================
    %% Stage 0
    %% =========================
    subgraph S0 [Stage 0 Define peptoid system]
        A0[Peptoid sequence]:::file
        A1[SMILES representation]:::file
        A2[Initial 3D structure generation]:::tool
        A0 --> A1 --> A2
        N0[Scope benchmarking]:::note
        N0 -.-> A0
    end

    %% =========================
    %% Stage 1
    %% =========================
    subgraph S1 [Stage 1 Geometry preparation]
        A2 --> B1[Generate initial coordinates]:::tool
        B1 --> B2[ORCA xTB geometry optimization]:::tool
        B2 --> B3[Optimized coordinates xyz]:::file
        B3 --> B4[PDB for atom identity]:::file
    end

    S1_OUT[Stage 1 output geometry]:::process
    B3 --> S1_OUT
    B4 --> S1_OUT

    %% =========================
    %% Stage 2
    %% =========================
    subgraph S2 [Stage 2 RESP charges]
        S1_OUT --> C1[NWChem ESP calculation]:::tool
        C1 --> C2[RESP charges]:::file
        N2[Cap groups neutral and fixed atom mapping]:::note
        N2 -.-> C2
    end

    S2_OUT[Stage 2 output charges]:::process
    C2 --> S2_OUT

    %% =========================
    %% Stage 3
    %% =========================
    subgraph S3 [Stage 3 Initial topology]
        S1_OUT --> D1[Antechamber residue prep]:::tool
        S2_OUT --> D1
        D1 --> D2[Prepi file]:::file
        D2 --> D3[tLeap topology build]:::tool
        D3 --> D4[Prmtop]:::file
        D3 --> D5[Inpcrd]:::file
        N3[Baseline parameters from GAFF2]:::note
        N3 -.-> D4
    end

    S3_OUT[Stage 3 topology]:::process
    D4 --> S3_OUT
    D5 --> S3_OUT

    %% =========================
    %% Stage 4
    %% =========================
    subgraph S4 [Stage 4 QM training set]
        S3_OUT --> E1[mdgx systematic sampling]:::tool
        E1 --> E2[QM input files]:::file
        E2 --> E3[QM single point calculations]:::tool
        E3 --> E4[QM outputs]:::file
        E4 --> E5[Training set energies and coordinates]:::file
        N4[QM engine swappable if parsing consistent]:::note
        N4 -.-> E3
    end

    S4_OUT[Stage 4 training set]:::process
    E5 --> S4_OUT

    %% =========================
    %% Stage 5
    %% =========================
    subgraph S5 [Stage 5 Parameter fitting]
        S3_OUT --> F1[mdgx fit dihedral terms]:::tool
        S4_OUT --> F1
        F1 --> F2[STEPs like frcmod]:::file
    end

    S5_OUT[New force field]:::process
    F2 --> S5_OUT

    %% =========================
    %% Stage 6
    %% =========================
    subgraph S6 [Stage 6 Benchmarking comparison]
        S5_OUT --> G1[Force field option STEPs like]:::process
        S3_OUT --> G2[Force field option GAFF2 baseline]:::process
        G3[Force field option OpenFF]:::process
        G4[Force field option CGenFF]:::process

        G1 --> H1[Run MD with OpenMM]:::tool
        G2 --> H1
        G3 --> H1
        G4 --> H1

        H1 --> I1[Trajectory set]:::file
        I1 --> I2[Common analysis metrics]:::tool
        I2 --> I3[Benchmark results]:::file

        N6[Compare relative behavior not accuracy]:::note
        N6 -.-> I3
    end

    S0 --> S1 --> S1_OUT --> S2 --> S2_OUT --> S3 --> S3_OUT --> S4 --> S4_OUT --> S5 --> S5_OUT --> S6
