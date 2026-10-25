"""Isolated CPU inference: upstream model, no training/telemetry dependencies."""
import json
from pathlib import Path
import sys
from types import SimpleNamespace

class Config(SimpleNamespace):
    def __contains__(self, key):
        return hasattr(self, key)


def main():
    vendor = Path(sys.argv[1]).resolve()
    sys.path.insert(0, str(vendor))
    import torch
    from transfer_model import TransferModel
    torch.set_num_threads(2)
    torch.manual_seed(0)
    payload = json.load(sys.stdin)
    config = Config(platform=Config(thermompnn_dir=str(vendor)), model=Config(
        hidden_dims=[64, 32], subtract_mut=True, num_final_layers=2, freeze_weights=True,
        load_pretrained=True, lightattn=True))
    model = TransferModel(config)
    checkpoint = torch.load(vendor / "models" / "thermoMPNN_default.pt", map_location="cpu", weights_only=True)
    state = {key[len("model."):]: value for key, value in checkpoint["state_dict"].items() if key.startswith("model.")}
    model.load_state_dict(state, strict=True)
    model.eval()
    sequence = payload["sequence"]
    pdb = {"name": "atlas_parent", "num_of_chains": 1, "seq": sequence, "seq_chain_A": sequence,
           "coords_chain_A": {f"{atom}_chain_A": coords for atom, coords in payload["coordinates"].items()}}
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    mutations = [SimpleNamespace(position=i, wildtype=wt, mutation=mut) for i, wt in enumerate(sequence) for mut in alphabet]
    with torch.inference_mode():
        predicted, _ = model([pdb], mutations)
    values = [float(row["ddG"].item()) for row in predicted]
    rows = [{mut: values[i*20+j] for j, mut in enumerate(alphabet)} for i in range(len(sequence))]
    json.dump({"sequence": sequence, "predictions": rows, "torch_version": torch.__version__}, sys.stdout, allow_nan=False)

if __name__ == "__main__":
    main()
