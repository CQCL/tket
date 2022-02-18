conan create --profile=tket recipes/tket
cd pytket
pip install -e .
cd ../
python kyriakos.py