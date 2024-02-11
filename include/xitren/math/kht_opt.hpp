#pragma once

#include <array>
#include <cstdint>
#include <cstring>
#include <iostream>

namespace xitren::math {

template <std::size_t Width, std::size_t Height>
union vault {
public:
    constexpr vault() : parts_{{}, {}}
    {
    }

    constexpr explicit vault(std::array<std::uint8_t, Width * Height>& image)
        : parts_{{}, image}
    {
    }

    constexpr explicit vault(std::array<std::uint8_t, Width * Height>&& image)
        : parts_{{}, std::move(image)}
    {
    }

    constexpr std::array<std::uint8_t, Width * Height>&
    image()
    {
        return parts_.image;
    }

    constexpr std::array<std::uint8_t, Width * Height>&
    mirror()
    {
        return parts_.mirror;
    }

    void
    print() {
        for (std::size_t i{1};i <= Width * Height * 2;i++) {
            std::cout << static_cast<int>(
                reinterpret_cast<std::uint8_t*>(all_.data())[i-1]) << "\t";
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
    struct tag_vault {
        std::array<std::uint8_t, Width * Height> mirror;
        std::array<std::uint8_t, Width * Height> image;
    } parts_;
};

template <std::size_t Width, std::size_t Height, typename Batch = std::uint32_t, bool Debug = false>
class kht {
    union image {
        std::array<std::uint8_t, Width * Height> all;
        std::array<std::array<Batch, Width / sizeof(Batch)>, Height> lines;
    };

    template <typename Type>
    union converter {
        std::array<std::uint8_t, sizeof(Type)> bytes;
        Type data;
    };

public:
    static constexpr void
    convert(vault<Width, Height>& input)
    {
        kht<Width, Height, Batch, Debug> r{input};
    }

    static constexpr void
    convert(std::array<std::uint8_t, Width * Height>& image,
            std::array<std::uint8_t, Width * Height>& mirror)
    {
        kht<Width, Height, Batch, Debug> r{image, mirror};
    }

protected:
    inline virtual void
    copy_memory(void* dst, const void* src, std::size_t size)
    {
        memcpy(dst, src, size);
    }

private:
    image* data_;
    image* mirror_;

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

    constexpr explicit kht(std::array<std::uint8_t, Width * Height>& image,
                           std::array<std::uint8_t, Width * Height>& image_mirror)
        : data_{reinterpret_cast<union image*>(&image)},
          mirror_{reinterpret_cast<union image*>(&image_mirror)}
    {
        if constexpr (Height <= 2) {
            mirror_first();
        }
        straight();
        if constexpr (Height > 2) {
            mirror();
        }
    }

    constexpr inline void straight()
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
        for (std::uint32_t i{0};i < half;i++) {
            auto& current{buffer[i]};

            Batch* current_row0{current[0].begin()};
            Batch* current_row1{current[1].begin()};
            Batch* low_row0{data_->lines[i].begin()};
            Batch* high_row0{reinterpret_cast<Batch *>(
                reinterpret_cast<std::uint8_t*>(data_->lines[i + half].begin()) + i)};
            Batch* high_row1{reinterpret_cast<Batch *>(
                reinterpret_cast<std::uint8_t*>(data_->lines[i + half].begin()) + i + 1)};

            const auto columns = (Width - (i + 1)) / sizeof(Batch);
            for (std::size_t j{0};j < columns;j++) {
                *(current_row0++) = *(low_row0) + *(high_row0++);
                *(current_row1++) = *(low_row0++) + *(high_row1++);
            }
            const auto columns_last = (Width / sizeof(Batch)) - columns;
            for (std::size_t j{0};j < columns_last;j++) {
                *(current_row0++) = *(low_row0);
                *(current_row1++) = *(low_row0++);
            }
        }
        copy_memory(data_->lines.data(), buffer.data(), size);
        if constexpr (Debug) {
            std::cout << "Straight Part Height: " << Height << std::endl;
            for (std::size_t i{1};i <= size;i++) {
                std::cout << static_cast<int>(
                    reinterpret_cast<std::uint8_t*>(buffer.data())[i-1]) << "\t";
                if (!(i % Width)) {
                    std::cout << std::endl;
                }
            }
            std::cout << std::endl;
        }
    }

    constexpr inline void mirror()
    {
        constexpr std::uint32_t half{Height / 2};
        constexpr std::uint32_t size{Height * Width};
        std::array<std::array<std::array<Batch, Width / sizeof(Batch)>, 2>, half> buffer;
        for (std::uint32_t i{0};i < half;i++) {
            auto& current{buffer[i]};

            Batch* current_row0{current[0].begin()};
            Batch* current_row1{current[1].begin()};
            Batch* low_row0{mirror_->lines[i].begin()};
            Batch* high_row0{reinterpret_cast<Batch *>(
                reinterpret_cast<std::uint8_t*>(mirror_->lines[i + half].begin()) + i)};
            Batch* high_row1{reinterpret_cast<Batch *>(
                reinterpret_cast<std::uint8_t*>(mirror_->lines[i + half].begin()) + i + 1)};

            const auto columns = (Width - (i + 1)) / sizeof(Batch);
            for (std::size_t j{0};j < columns;j++) {
                *(current_row0++) = *(low_row0) + *(high_row0++);
                *(current_row1++) = *(low_row0++) + *(high_row1++);
            }
            const auto columns_last = (Width / sizeof(Batch)) - columns;
            for (std::size_t j{0};j < columns_last;j++) {
                *(current_row0++) = *(low_row0);
                *(current_row1++) = *(low_row0++);
            }
        }
        copy_memory(mirror_->lines.data(), buffer.data(), size);
        if constexpr (Debug) {
            std::cout << "Mirror Part Height: " << Height << std::endl;
            for (std::size_t i{1};i <= size;i++) {
                std::cout << static_cast<int>(
                    reinterpret_cast<std::uint8_t*>(buffer.data())[i-1]) << "\t";
                if (!(i % Width)) {
                    std::cout << std::endl;
                }
            }
            std::cout << std::endl;
        }
    }

    constexpr inline void mirror_first()
    {
        constexpr std::uint32_t half{Height / 2};
        constexpr std::uint32_t size{Height * Width};
        std::array<std::array<std::array<Batch, Width / sizeof(Batch)>, 2>, half> buffer;
        for (std::uint32_t i{0};i < half;i++) {
            auto& current{buffer[i]};
            auto k = ((half - (1 + i)) << 1);

            Batch* current_row0{current[0].begin()};
            Batch* current_row1{current[1].begin()};
            Batch* low_row0{data_->lines[k + 1].begin()};
            Batch* high_row0{reinterpret_cast<Batch *>(
                reinterpret_cast<std::uint8_t*>(data_->lines[k].begin()) + i)};
            Batch* high_row1{reinterpret_cast<Batch *>(
                reinterpret_cast<std::uint8_t*>(data_->lines[k].begin()) + i + 1)};

            const auto columns = (Width - (i + 1)) / sizeof(Batch);
            for (std::size_t j{0};j < columns;j++) {
                *(current_row0++) = *(low_row0) + *(high_row0++);
                *(current_row1++) = *(low_row0++) + *(high_row1++);
            }
            const auto columns_last = (Width / sizeof(Batch)) - columns;
            for (std::size_t j{0};j < columns_last;j++) {
                *(current_row0++) = *(low_row0);
                *(current_row1++) = *(low_row0++);
            }
        }
        copy_memory(mirror_->lines.data(), buffer.data(), size);
        if constexpr (Debug) {
            std::cout << "Mirror Part Height: " << Height << std::endl;
            for (std::size_t i{1};i <= size;i++) {
                std::cout << static_cast<int>(
                    reinterpret_cast<std::uint8_t*>(buffer.data())[i-1]) << "\t";
                if (!(i % Width)) {
                    std::cout << std::endl;
                }
            }
            std::cout << std::endl;
        }
    }
};

}
