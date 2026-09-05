#include "raym0nade/gpu/vulkan_capabilities.hpp"
#include "raym0nade/gpu/vulkan_self_test.hpp"

#include <cstdint>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

double gibibytes(std::uint64_t bytes) noexcept {
    constexpr double bytesPerGibibyte = 1024.0 * 1024.0 * 1024.0;
    return static_cast<double>(bytes) / bytesPerGibibyte;
}

void printVersion(const raym0nade::gpu::VulkanVersion& version) {
    std::cout << version.major << '.' << version.minor << '.' << version.patch;
}

struct Options {
    bool runSelfTest{false};
    bool requestValidation{false};
};

void printUsage() {
    std::cout << "Usage: raym0nade_gpu_probe [--self-test] [--validation]\n"
                 "  --self-test   Build a triangle BLAS/TLAS and execute deterministic Ray Queries.\n"
                 "  --validation  Enable the Khronos validation layer; implies --self-test.\n";
}

Options parseOptions(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument{argv[index]};
        if (argument == "--self-test") {
            options.runSelfTest = true;
        } else if (argument == "--validation") {
            options.runSelfTest = true;
            options.requestValidation = true;
        } else if (argument == "--help" || argument == "-h") {
            printUsage();
            std::exit(0);
        } else {
            throw std::invalid_argument("Unknown option: " + argument);
        }
    }
    return options;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parseOptions(argc, argv);
        const raym0nade::gpu::VulkanVersion loader = raym0nade::gpu::vulkanLoaderVersion();
        std::cout << "Vulkan loader: ";
        printVersion(loader);
        std::cout << '\n';

        const std::vector<raym0nade::gpu::VulkanDeviceCapabilities> devices =
            raym0nade::gpu::enumerateVulkanDevices();
        if (devices.empty()) {
            std::cerr << "No Vulkan physical devices were found.\n";
            return options.runSelfTest ? 77 : 2;
        }

        bool foundSupportedAmdDevice = false;
        for (std::size_t index = 0; index < devices.size(); ++index) {
            const auto& device = devices[index];
            std::cout << "\nDevice " << index << ": " << device.deviceName << '\n'
                      << "  Vendor/device: 0x" << std::hex << std::setw(4) << std::setfill('0')
                      << device.vendorId << "/0x" << std::setw(4) << device.deviceId << std::dec
                      << std::setfill(' ') << '\n'
                      << "  API: ";
            printVersion(device.apiVersion);
            std::cout << '\n'
                      << "  Driver: " << device.driverName << " (" << device.driverInfo << ")\n"
                      << "  Type: " << (device.integrated ? "integrated" : "discrete-or-other") << '\n'
                      << "  Device-local heaps: " << std::fixed << std::setprecision(2)
                      << gibibytes(device.deviceLocalMemoryBytes) << " GiB\n"
                      << "  Compute queue: " << (device.hasComputeQueue ? "yes" : "no");
            if (device.hasComputeQueue) {
                std::cout << " (family " << device.computeQueueFamily << ')';
            }
            std::cout << '\n'
                      << "  Subgroup size: " << device.subgroupSize << '\n'
                      << "  Buffer device address: " << (device.bufferDeviceAddress ? "yes" : "no") << '\n'
                      << "  Acceleration structure: " << (device.accelerationStructure ? "yes" : "no") << '\n'
                      << "  Ray Query: " << (device.rayQuery ? "yes" : "no") << '\n'
                      << "  Ray Tracing Pipeline: " << (device.rayTracingPipeline ? "yes" : "no") << '\n'
                      << "  Maximum AS primitives: "
                      << device.maximumAccelerationStructurePrimitiveCount << '\n';

            if (device.supportsRayQueryBackend()) {
                foundSupportedAmdDevice = true;
                std::cout << "  Raym0nade AMD Ray Query backend: supported\n";
            } else {
                std::cout << "  Raym0nade AMD Ray Query backend: unsupported\n";
                for (const std::string& requirement : device.missingRequirements) {
                    std::cout << "    Missing: " << requirement << '\n';
                }
                if (!device.isAmd()) {
                    std::cout << "    Missing: AMD vendor device\n";
                }
            }
        }

        if (!foundSupportedAmdDevice) {
            return options.runSelfTest ? 77 : 3;
        }
        if (!options.runSelfTest) {
            return 0;
        }

        const auto selfTest = raym0nade::gpu::runVulkanRayQuerySelfTest(
            std::filesystem::path{RAYM0NADE_RAY_QUERY_SELF_TEST_SHADER},
            options.requestValidation);
        std::cout << '\n' << raym0nade::gpu::formatVulkanRayQuerySelfTestResult(selfTest);
        if (options.requestValidation && !selfTest.validationEnabled) {
            std::cerr << "The requested Khronos validation layer was unavailable.\n";
            return 5;
        }
        return selfTest.passed ? 0 : 4;
    } catch (const std::exception& error) {
        std::cerr << "GPU capability probe failed: " << error.what() << '\n';
        return 1;
    }
}
