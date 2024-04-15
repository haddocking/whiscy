import shutil
from pathlib import Path

import pytest


@pytest.fixture
def scratch_path(request):
    basename = request.param

    path = Path(Path(__file__).parent, basename)
    path.mkdir(parents=True, exist_ok=True)

    yield path

    shutil.rmtree(path)
