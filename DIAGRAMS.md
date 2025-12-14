# ONT-SMA-seq Diagrams

Quick reference diagrams for the ONT-SMA-seq pipeline.

## Pipeline Architecture

```mermaid
graph TB
    subgraph "Step 1: Database Initialization"
        A1[mkdb.py] -->|Creates| A2[SMA_exp_id.db]
        A2 -->|Contains| A3[Schema: Reads, Mods, Exp, Refseq]
    end
    
    subgraph "Step 2: Input Standardization"
        B1[inputInit.py] -->|Symlinks| B2[Input/ Directory]
        B2 -->|Contains| B3[Standardized File Paths]
    end
    
    subgraph "Step 3: Data Processing"
        C1[ingest.py] -->|Processes| C2[Read Data]
        C2 -->|Populates| C3[Database Tables]
        C2 -->|Creates| C4[Tagged BAM Output]
    end
    
    A2 --> C1
    B2 --> C1
    
    style A1 fill:#0077b6,stroke:#004080,stroke-width:2px,color:#fff
    style B1 fill:#0077b6,stroke:#004080,stroke-width:2px,color:#fff
    style C1 fill:#0077b6,stroke:#004080,stroke-width:2px,color:#fff
    style A2 fill:#fd7e14,stroke:#dc6308,stroke-width:2px,color:#fff
    style B2 fill:#fd7e14,stroke:#dc6308,stroke-width:2px,color:#fff
    style C3 fill:#28a745,stroke:#1e7e34,stroke-width:2px,color:#fff
    style C4 fill:#28a745,stroke:#1e7e34,stroke-width:2px,color:#fff
```

## Input/Output Summary

```mermaid
flowchart LR
    subgraph Inputs
        I1[Raw BAM<br/>filename with metadata]
        I2[Pod5 Directory<br/>end reason data]
        I3[Reference FASTA<br/>2 sequences]
        I4[Experiment ID]
    end
    
    subgraph Pipeline
        P1[mkdb.py]
        P2[inputInit.py]
        P3[ingest.py]
    end
    
    subgraph Outputs
        O1[SQLite Database<br/>Complete data]
        O2[Tagged BAM<br/>with ER field]
        O3[Standardized Input/<br/>directory structure]
    end
    
    I4 --> P1 --> O1
    I1 --> P2
    I2 --> P2
    I3 --> P2
    P2 --> O3
    O3 --> P3
    O1 --> P3
    P3 --> O1
    P3 --> O2
    
    style Inputs fill:#fd7e14,stroke:#dc6308,stroke-width:2px,color:#fff
    style Pipeline fill:#0077b6,stroke:#004080,stroke-width:2px,color:#fff
    style Outputs fill:#28a745,stroke:#1e7e34,stroke-width:2px,color:#fff
```

## Data Structure: Database Schema

```mermaid
erDiagram
    Reads {
        TEXT uniq_id PK
        TEXT exp_id FK
        TEXT refseq_id FK
        TEXT read_id
        TEXT readseq
        INTEGER readlen
        TEXT model_tier
        TEXT model_ver
        INTEGER trim
        INTEGER mod_bitflag FK
        INTEGER ed
        REAL q_bc
        REAL q_ld
        TEXT ER
    }
    
    Mods {
        INTEGER mod_bitflag PK
        TEXT mods
    }
    
    Exp {
        TEXT exp_id PK
        TEXT exp_desc
    }
    
    Refseq {
        TEXT refseq_id PK
        TEXT refseq
        INTEGER reflen
    }
    
    Reads }|--|| Exp : "exp_id"
    Reads }o--|| Refseq : "refseq_id"
    Reads }|--|| Mods : "mod_bitflag"
```

## Read Processing Flow

```mermaid
flowchart TD
    Start([Read from BAM]) --> Parse[Parse Read Data]
    Parse --> LookupER[Lookup End Reason<br/>from Pod5]
    LookupER --> TagER[Add ER Tag]
    TagER --> WriteOut[Write to Output BAM]
    WriteOut --> CheckLen{Read Length<br/>Matches Reference?}
    
    CheckLen -->|Yes| CalcED[Calculate Edit Distance<br/>Levenshtein]
    CheckLen -->|No| SetNull[ed = NULL<br/>q_ld = NULL<br/>refseq_id = NULL]
    
    CalcED --> CalcQLD[Calculate q_ld]
    CalcQLD --> SetRef[Set refseq_id]
    
    SetRef --> CalcQBC[Calculate q_bc<br/>Basecall Quality]
    SetNull --> CalcQBC
    
    CalcQBC --> GenID[Generate uniq_id]
    GenID --> InsertDB[(Insert into<br/>Reads Table)]
    InsertDB --> Next{More Reads?}
    
    Next -->|Yes| Start
    Next -->|No| Done([Complete])
    
    style Start fill:#0077b6,stroke:#004080,stroke-width:2px,color:#fff
    style Done fill:#28a745,stroke:#1e7e34,stroke-width:2px,color:#fff
    style InsertDB fill:#fd7e14,stroke:#dc6308,stroke-width:2px,color:#fff
```

## Data Transformations

### Input BAM Filename → Metadata

```
Filename: EXP001_h_v5.2.0_1_6mA.bam
         ↓
Parsed Metadata:
  • exp_id: EXP001
  • model_tier: h (high accuracy)
  • model_ver: 5.2.0
  • trim: 1 (trimmed)
  • mods_str: 6mA
  • mod_bitflag: 1
```

### Read → Database Record

```
Input Read:
  • read_id: abc123...
  • sequence: ATCG...
  • quality: [30, 32, 28...]
  • length: 1500

Processing:
  • Match to reference (length 1450 ± 150) ✓
  • Calculate ed = 15 (Levenshtein)
  • Calculate q_bc = 31.5
  • Calculate q_ld = 19.2
  • Generate uniq_id: EXP001h520t1m1_a3f2d8e1

Database Record:
  uniq_id: EXP001h520t1m1_a3f2d8e1
  exp_id: EXP001
  refseq_id: reference_short
  read_id: abc123...
  readseq: ATCG...
  readlen: 1500
  model_tier: h
  model_ver: 5.2.0
  trim: 1
  mod_bitflag: 1
  ed: 15
  q_bc: 31.5
  q_ld: 19.2
  ER: signal_positive
```

## Modification Bitflags

```mermaid
graph LR
    subgraph "Base Modifications"
        M0[non = 0]
        M1[6mA = 1]
        M2[5mCG_5hmCG = 2]
        M3[5mC_5hmC = 4]
        M4[4mC_5mC = 8]
        M5[5mC = 16]
    end
    
    subgraph "Combinations 6mA + C-mod"
        C1[6mA + 5mCG_5hmCG = 3]
        C2[6mA + 5mC_5hmC = 5]
        C3[6mA + 4mC_5mC = 9]
        C4[6mA + 5mC = 17]
    end
    
    M1 --> C1
    M2 --> C1
    M1 --> C2
    M3 --> C2
    M1 --> C3
    M4 --> C3
    M1 --> C4
    M5 --> C4
    
    style M0 fill:#6c757d,stroke:#495057,stroke-width:2px,color:#fff
    style M1 fill:#0077b6,stroke:#004080,stroke-width:2px,color:#fff
    style M2 fill:#ffc107,stroke:#d39e00,stroke-width:2px,color:#000
    style M3 fill:#ffc107,stroke:#d39e00,stroke-width:2px,color:#000
    style M4 fill:#ffc107,stroke:#d39e00,stroke-width:2px,color:#000
    style M5 fill:#ffc107,stroke:#d39e00,stroke-width:2px,color:#000
    style C1 fill:#28a745,stroke:#1e7e34,stroke-width:2px,color:#fff
    style C2 fill:#28a745,stroke:#1e7e34,stroke-width:2px,color:#fff
    style C3 fill:#28a745,stroke:#1e7e34,stroke-width:2px,color:#fff
    style C4 fill:#28a745,stroke:#1e7e34,stroke-width:2px,color:#fff
```

## Directory Structure Before/After

```mermaid
graph TB
    subgraph "Before Pipeline"
        B1[project/]
        B2[data/raw.bam]
        B3[data/pod5/]
        B4[data/ref.fa]
        B1 --> B2
        B1 --> B3
        B1 --> B4
    end
    
    subgraph "After Pipeline"
        A1[project/]
        A2[SMA_exp.db]
        A3[Input/exp*.bam symlink]
        A4[Input/exp_pod5/ symlink]
        A5[Input/exp.fa symlink]
        A6[Output/exp.bam tagged]
        A1 --> A2
        A1 --> A3
        A1 --> A4
        A1 --> A5
        A1 --> A6
    end
    
    style A2 fill:#fd7e14,stroke:#dc6308,stroke-width:2px,color:#fff
    style A6 fill:#28a745,stroke:#1e7e34,stroke-width:2px,color:#fff
```

## Metric Calculations

### Basecall Quality (q_bc)

```
Input: Quality scores [Q1, Q2, Q3, ..., Qn]

Step 1: Convert to error probabilities
  P_error = 10^(-Q/10)

Step 2: Calculate mean error
  mean_P = Σ(P_error) / n

Step 3: Convert back to quality
  q_bc = -10 * log10(mean_P)
```

### Levenshtein Quality (q_ld)

```
Input: 
  - ed (edit distance)
  - L (reference length)

Step 1: Calculate ratio
  ratio = ed / L

Step 2: Apply bounds
  bounded = min(max(1/L², ratio), 1)

Step 3: Convert to quality
  q_ld = -10 * log10(bounded)
```

### Unique ID Generation

```
Components:
  - exp_id: EXP001
  - tier: h
  - ver: 5.2.0 → 520
  - trim: 1
  - mod_bitflag: 1
  - read_hash: MD5(read_id)[:8]

Format: {exp_id}{tier}{ver}t{trim}m{mod}_{hash}
Result: EXP001h520t1m1_a3f2d8e1
```

## Reference Matching Logic

```mermaid
flowchart TD
    Start([Read Length L]) --> Check1{L in long range?<br/>Llong ± 150}
    Check1 -->|Yes| Match1[Match: Long Reference]
    Check1 -->|No| Check2{L in short range?<br/>Lshort ± 150}
    Check2 -->|Yes| Match2[Match: Short Reference]
    Check2 -->|No| NoMatch[No Match<br/>refseq_id = NULL]
    
    Match1 --> Calc[Calculate Metrics:<br/>ed, q_ld]
    Match2 --> Calc
    NoMatch --> Skip[Skip Metrics:<br/>ed = NULL<br/>q_ld = NULL]
    
    Calc --> DB[(Store in Database)]
    Skip --> DB
    
    style Start fill:#0077b6,stroke:#004080,stroke-width:2px,color:#fff
    style Match1 fill:#28a745,stroke:#1e7e34,stroke-width:2px,color:#fff
    style Match2 fill:#28a745,stroke:#1e7e34,stroke-width:2px,color:#fff
    style NoMatch fill:#dc3545,stroke:#bd2130,stroke-width:2px,color:#fff
    style DB fill:#fd7e14,stroke:#dc6308,stroke-width:2px,color:#fff
```
