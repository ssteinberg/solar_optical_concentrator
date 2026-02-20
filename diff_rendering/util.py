import csv
import json
import time
from pathlib import Path
from typing import Any

import mitsuba as mi
from matplotlib.figure import Figure

OUTPUT_FOLDER_NAME = 'output'

def create_and_get_output_path(project_name: str) -> Path:
    """
    Creates a folder to hold the output of the project.
    :param project_name: The name of the project.
    :return: The path to the created, output folder.
    """
    output_path = Path(f'{OUTPUT_FOLDER_NAME}/{project_name}/{time.strftime("%Y-%m-%d_%H-%M-%S")}')
    output_path.mkdir(parents=True, exist_ok=True)
    return output_path

def write_data_json(data_dict: dict, output_path: Path, file_name: str) -> None:
    """
    Write the provided dictionary to a JSON file.
    :param data_dict: The data to write.
    :param output_path: The output folder for the run.
    :param file_name: The JSON file name.
    """
    with open(output_path.joinpath(f'{file_name}.json'), 'w') as f:
        json.dump(data_dict, f, indent=4)

def write_render(render: Any, output_path: Path, file_name: str) -> None:
    """
    Write the provided mitsuba render to a PNG file.
    :param render: The render to write.
    :param output_path: The output folder for the run.
    :param file_name: The PNG file name.
    """
    mi.util.write_bitmap(str(output_path.joinpath(f'{file_name}.png')), render)

def write_figure(fig: Figure, output_path: Path, file_name: str) -> None:
    """
    Write the provided matplotlib figure to a PNG file.
    :param fig: The figure to write.
    :param output_path: The output folder for the run.
    :param file_name: The PNG file name.
    """
    fig.savefig(output_path.joinpath(f'{file_name}.png'), bbox_inches="tight", pad_inches=0, dpi=1000)

def write_mirror(x: Any, y: Any, output_path: Path, file_name: str) -> None:
    """
    Write the provided mirror vertices to a CSV file.
    :param x: The x coordinates of the vertices.
    :param y: The y coordinates of the vertices.
    :param output_path: The output folder for the run.
    :param file_name: The JSON file name.
    """
    with open(output_path.joinpath(f'{file_name}.csv'), 'w', newline='') as f:
        csv.writer(f).writerows(zip(x, y))
