# Quadratube

## Introduction
Disloactions pair on quadrangular nanotube are studied here.

![triangular grid](doc/fig/grid_3.svg)

## Model Details
See documentation [here](doc/md/model.md).

## Code Specification
See documentation [here](doc/md/code.md).

## Brief Introduction Of Kokkos
See documentation [here](doc/md/kokkos.md).

## Ovito Column Mapping Preset
To load ovito column mapping preset, add [this](etc/ovito_preset.txt) into your Ovito.conf file. This file always locates at `~/.config/Ovito/Ovito.conf` on Linux. Or you can manually set like this:
```json
{
    "g_curv": "Gaussian Curvature",
    "m_curv": "Mean Curvature"
}
```

## Test
`test-inl.h` is for code test. Write `main` file like this

```cpp
#include "test-inl.h"

int main() {
  test_derivative(); // test function
}
```

## Compile
If you want to use model 3 and want higher speed, run
```sh
cmake .. -DModel3=ON -DebugType=OFF
```

## Bugs
See documentation [here](doc/md/bugs.md).