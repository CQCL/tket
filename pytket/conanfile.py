from conan import ConanFile
from conan.tools.cmake import CMakeToolchain, CMake, cmake_layout, CMakeDeps


class pytketRecipe(ConanFile):
    name = "pytket"
    version = "1.0.0"
    package_type = "application"

    # Binary configuration
    settings = "os", "compiler", "build_type", "arch"

    # Sources are located in the same place as this recipe, copy them to the recipe
    exports_sources = "CMakeLists.txt", "binders/*"

    def layout(self):
        cmake_layout(self)

    def generate(self):
        deps = CMakeDeps(self)
        deps.generate()
        tc = CMakeToolchain(self)
        tc.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        cmake = CMake(self)
        cmake.install()

    def requirements(self):
        self.requires("tket/1.3.0@tket/stable")
        self.requires("tklog/0.3.3@tket/stable")
        self.requires("tkrng/0.3.3@tket/stable")
        self.requires("tkassert/0.3.4@tket/stable")
        self.requires("tkwsm/0.3.8@tket/stable")
        self.requires("tktokenswap/0.3.8@tket/stable")
        self.requires("symengine/0.11.2")
        self.requires("gmp/6.3.0")
        self.requires("pybind11/2.12.0")
        self.requires("nlohmann_json/3.11.3")
        self.requires("pybind11_json/0.2.14")
