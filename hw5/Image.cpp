// define our communication data type

struct IMAGE
{
    unsigned char image[SIZE];

    IMAGE(void)
    {
        for (int i = 0; i < SIZE; i++)
            image[i] = 0;
    }

    IMAGE& operator=(const IMAGE& copy)
    {
        for (int i = 0; i < SIZE; i++)
            image[i] = copy.image[i];
        return *this;
    }

    unsigned char& operator[](const int index)
    {
        return image[index];
    }

    operator unsigned char*()
    {
        return image;
    }
};
