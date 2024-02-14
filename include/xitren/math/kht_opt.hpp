#pragma once

#include <array>
#include <cstdint>
#include <cstring>
#include <iostream>

namespace xitren::math {

template <std::size_t Width, std::size_t Height>
union vault {
public:
    /**
     * Constructs a new instance of the vault class.
     * The parts array is initialized to all zeros.
     */
    constexpr vault() : parts_{{}, {}} {}

    /**
     * Constructs a new instance of the vault class.
     * The parts array is initialized to the contents of the specified image array.
     * @param image The array containing the image data.
     */
    constexpr explicit vault(std::array<std::uint8_t, Width * Height>& image) : parts_{{}, image} {}

    /**
     * Constructs a new instance of the vault class.
     * The parts array is initialized to the contents of the specified rvalue image array.
     * @param image The rvalue array containing the image data.
     */
    constexpr explicit vault(std::array<std::uint8_t, Width * Height>&& image) : parts_{{}, std::move(image)} {}

    /**
     * Returns a reference to the image part of the vault.
     * @return A reference to the image part of the vault.
     */
    constexpr std::array<std::uint8_t, Width * Height>&
    image()
    {
        return parts_.image;
    }

    /**
     * Returns a reference to the mirror part of the vault.
     * @return A reference to the mirror part of the vault.
     */
    constexpr std::array<std::uint8_t, Width * Height>&
    mirror()
    {
        return parts_.mirror;
    }

    /**
     * Prints the contents of the vault to the console.
     */
    void
    print()
    {
        for (std::size_t i{1}; i <= Width * Height * 2; i++) {
            std::cout << static_cast<int>(reinterpret_cast<std::uint8_t*>(all_.data())[i - 1]) << "\t";
            if (!(i % Width)) {
                std::cout << std::endl;
            }
            if (i == (Width * Height))
                std::cout << std::endl;
        }
        std::cout << std::endl;
    }

private:
    std::array<std::uint8_t, Width * Height * 2> all_;

    /**
     * The parts of the vault.
     * The first element of the array is the mirror part, and the second element is the image part.
     */
    struct tag_vault {
        std::array<std::uint8_t, Width * Height> mirror;
        std::array<std::uint8_t, Width * Height> image;
    } parts_;
};

template <std::size_t Width, std::size_t Height, typename Batch = std::uint32_t, bool Debug = false>
class kht {
    union image {
        std::array<std::uint8_t, Width * Height>                     all;
        std::array<std::array<Batch, Width / sizeof(Batch)>, Height> lines;
    };

    /**
     * A union that represents a converter between a data type and an array of bytes.
     * The converter can be used to convert between the data type and an array of bytes, or between two arrays of bytes.
     * The converter stores the data in the form of an array of bytes, and provides a way to access the data as a data
     * type. The converter is implemented as a union, which allows the data to be stored in a single memory location,
     * and accessed as either a data type or an array of bytes.
     * @tparam Type The data type to be converted.
     */
    template <typename Type>
    union converter {
        std::array<std::uint8_t, sizeof(Type)> bytes;
        /**
         * The data of the converter.
         */
        Type data;
    };

public:
    /**
     * A static method that converts the input data in the vault class.
     * @param input the input data in the vault class.
     */
    static constexpr void
    convert(vault<Width, Height>& input)
    {
        kht<Width, Height, Batch, Debug> r{input};
    }

    /**
     * A static method that converts the input data in the std::array class.
     * @param image the input data in the std::array class.
     * @param mirror the mirror input data in the std::array class.
     */
    static constexpr void
    convert(std::array<std::uint8_t, Width * Height>& image, std::array<std::uint8_t, Width * Height>& mirror)
    {
        kht<Width, Height, Batch, Debug> r{image, mirror};
    }

protected:
    /**
     * An inline virtual function that copies the memory from one location to another.
     * @param dst the destination location.
     * @param src the source location.
     * @param size the size of the memory to be copied.
     */
    inline virtual void
    copy_memory(void* dst, void const* src, std::size_t size)
    {
        memcpy(dst, src, size);
    }

private:
    /**
     * A pointer to the data in the vault class.
     */
    image* data_;
    /**
     * A pointer to the mirror data in the vault class.
     */
    image* mirror_;

    /**
     * A constructor that takes the input data in the vault class as input.
     * @param input the input data in the vault class.
     */
    constexpr explicit kht(vault<Width, Height>& input)
        : data_{reinterpret_cast<union image*>(&input.image())},
          mirror_{reinterpret_cast<union image*>(&input.mirror())}
    {
        if constexpr (Height <= 2) {
            mirror_first();
        }
        straight();
        if constexpr (Height > 2) {
            mirror();
        }
    }

    /**
     * A constructor that takes the input data in the std::array class as input.
     * @param image the input data in the std::array class.
     * @param mirror the mirror input data in the std::array class.
     */
    constexpr explicit kht(std::array<std::uint8_t, Width * Height>& image,
                           std::array<std::uint8_t, Width * Height>& image_mirror)
        : data_{reinterpret_cast<union image*>(&image)}, mirror_{reinterpret_cast<union image*>(&image_mirror)}
    {
        if constexpr (Height <= 2) {
            mirror_first();
        }
        straight();
        if constexpr (Height > 2) {
            mirror();
        }
    }

    /**
     * A function that performs the straight part of the KHT algorithm.
     */
    inline constexpr void
    straight()
    {
        constexpr std::uint32_t half{Height / 2};
        constexpr std::uint32_t size{Height * Width};
        if constexpr (half > 1) {
            kht<Width, half, Batch, Debug>::convert(
                *reinterpret_cast<std::array<std::uint8_t, Width * half>*>(data_),
                *(reinterpret_cast<std::array<std::uint8_t, Width * half>*>(mirror_) + 1));
            kht<Width, half, Batch, Debug>::convert(
                *(reinterpret_cast<std::array<std::uint8_t, Width * half>*>(data_) + 1),
                *reinterpret_cast<std::array<std::uint8_t, Width * half>*>(mirror_));
        }
        std::array<std::array<std::array<Batch, Width / sizeof(Batch)>, 2>, half> buffer;
        for (std::uint32_t i{0}; i < half; i++) {
            auto& current{buffer[i]};

            Batch* current_row0{current[0].begin()};
            Batch* current_row1{current[1].begin()};
            Batch* low_row0{data_->lines[i].begin()};
            Batch* high_row0{
                reinterpret_cast<Batch*>(reinterpret_cast<std::uint8_t*>(data_->lines[i + half].begin()) + i)};
            Batch* high_row1{
                reinterpret_cast<Batch*>(reinterpret_cast<std::uint8_t*>(data_->lines[i + half].begin()) + i + 1)};

            auto const columns = (Width - (i + 1)) / sizeof(Batch);
            for (std::size_t j{0}; j < columns; j++) {
                *(current_row0++) = *(low_row0) + *(high_row0++);
                *(current_row1++) = *(low_row0++) + *(high_row1++);
            }
            auto const columns_last = (Width / sizeof(Batch)) - columns;
            for (std::size_t j{0}; j < columns_last; j++) {
                *(current_row0++) = *(low_row0);
                *(current_row1++) = *(low_row0++);
            }
        }
        copy_memory(data_->lines.data(), buffer.data(), size);
        if constexpr (Debug) {
            std::cout << "Straight Part Height: " << Height << std::endl;
            for (std::size_t i{1}; i <= size; i++) {
                std::cout << static_cast<int>(reinterpret_cast<std::uint8_t*>(buffer.data())[i - 1]) << "\t";
                if (!(i % Width)) {
                    std::cout << std::endl;
                }
            }
            std::cout << std::endl;
        }
    }

    /**
     * A function that performs the mirror part of the KHT algorithm.
     */
    inline constexpr void
    mirror()
    {
        constexpr std::uint32_t                                                   half{Height / 2};
        constexpr std::uint32_t                                                   size{Height * Width};
        std::array<std::array<std::array<Batch, Width / sizeof(Batch)>, 2>, half> buffer;
        for (std::uint32_t i{0}; i < half; i++) {
            auto& current{buffer[i]};

            Batch* current_row0{current[0].begin()};
            Batch* current_row1{current[1].begin()};
            Batch* low_row0{mirror_->lines[i].begin()};
            Batch* high_row0{
                reinterpret_cast<Batch*>(reinterpret_cast<std::uint8_t*>(mirror_->lines[i + half].begin()) + i)};
            Batch* high_row1{
                reinterpret_cast<Batch*>(reinterpret_cast<std::uint8_t*>(mirror_->lines[i + half].begin()) + i + 1)};

            auto const columns = (Width - (i + 1)) / sizeof(Batch);
            for (std::size_t j{0}; j < columns; j++) {
                *(current_row0++) = *(low_row0) + *(high_row0++);
                *(current_row1++) = *(low_row0++) + *(high_row1++);
            }
            auto const columns_last = (Width / sizeof(Batch)) - columns;
            for (std::size_t j{0}; j < columns_last; j++) {
                *(current_row0++) = *(low_row0);
                *(current_row1++) = *(low_row0++);
            }
        }
        copy_memory(mirror_->lines.data(), buffer.data(), size);
        if constexpr (Debug) {
            std::cout << "Mirror Part Height: " << Height << std::endl;
            for (std::size_t i{1}; i <= size; i++) {
                std::cout << static_cast<int>(reinterpret_cast<std::uint8_t*>(buffer.data())[i - 1]) << "\t";
                if (!(i % Width)) {
                    std::cout << std::endl;
                }
            }
            std::cout << std::endl;
        }
    }

    /**
     * A function that performs the mirror first part of the KHT algorithm.
     */
    inline constexpr void
    mirror_first()
    {
        constexpr std::uint32_t                                                   half{Height / 2};
        constexpr std::uint32_t                                                   size{Height * Width};
        std::array<std::array<std::array<Batch, Width / sizeof(Batch)>, 2>, half> buffer;
        for (std::uint32_t i{0}; i < half; i++) {
            auto& current{buffer[i]};
            auto  k = ((half - (1 + i)) << 1);

            Batch* current_row0{current[0].begin()};
            Batch* current_row1{current[1].begin()};
            Batch* low_row0{data_->lines[k + 1].begin()};
            Batch* high_row0{reinterpret_cast<Batch*>(reinterpret_cast<std::uint8_t*>(data_->lines[k].begin()) + i)};
            Batch* high_row1{
                reinterpret_cast<Batch*>(reinterpret_cast<std::uint8_t*>(data_->lines[k].begin()) + i + 1)};

            auto const columns = (Width - (i + 1)) / sizeof(Batch);
            for (std::size_t j{0}; j < columns; j++) {
                *(current_row0++) = *(low_row0) + *(high_row0++);
                *(current_row1++) = *(low_row0++) + *(high_row1++);
            }
            auto const columns_last = (Width / sizeof(Batch)) - columns;
            for (std::size_t j{0}; j < columns_last; j++) {
                *(current_row0++) = *(low_row0);
                *(current_row1++) = *(low_row0++);
            }
        }
        copy_memory(mirror_->lines.data(), buffer.data(), size);
        if constexpr (Debug) {
            std::cout << "Mirror Part Height: " << Height << std::endl;
            for (std::size_t i{1}; i <= size; i++) {
                std::cout << static_cast<int>(reinterpret_cast<std::uint8_t*>(buffer.data())[i - 1]) << "\t";
                if (!(i % Width)) {
                    std::cout << std::endl;
                }
            }
            std::cout << std::endl;
        }
    }
};

}    // namespace xitren::math
