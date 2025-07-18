from conan import ConanFile
from conan.tools.cmake import CMake, CMakeDeps, CMakeToolchain, cmake_layout


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
        self.requires("gmp/6.3.0")
        self.requires("nlohmann_json/3.12.0")
        self.requires("nanobind/tci-2.8.0@tket/stable")
        self.requires("symengine/tci-0.14.0.1@tket/stable")
        self.requires("tkassert/0.3.4@tket/stable")
        self.requires("tket/2.1.35@tket/stable")
        self.requires("tklog/0.3.3@tket/stable")
        self.requires("tkrng/0.3.3@tket/stable")
        self.requires("tktokenswap/0.3.11@tket/stable")
        self.requires("tkwsm/0.3.11@tket/stable")
