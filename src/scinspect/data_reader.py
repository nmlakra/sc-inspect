from scanpy import read_h5ad
from scanpy import AnnData


def read_data(filepath: str) -> AnnData:
    adata = read_h5ad(filepath, backed="r+")
    return adata
