#include <iostream>
#include <random>
#include <map>
#include <array>
#include <stdint.h>
#include <functional>
#include <algorithm>
#include <sstream>
#include <numeric>
#include <ostream>
#include <fstream>

#define VERBOSE 0

std::random_device rd;  // a seed source for the random number engine
std::mt19937 gen(rd()); // mersenne_twister_engine seeded with rd()
std::map<int, std::map<int, std::uniform_int_distribution<>>> randomDistributions;

int rollDice(int min, int max)
{
    std::uniform_int_distribution<> *rng = nullptr;
    auto minMap = randomDistributions.find(min);
    if (minMap != randomDistributions.end())
    {
        auto maxMap = minMap->second.find(max);
        if (maxMap != minMap->second.end())
        {
#if VERBOSE
            std::cout << "Reusing already existing random distribution: [" << min << "," << max << "]" << std::endl;
#endif
            rng = &(maxMap->second);
        }
        else
        {
#if VERBOSE
            std::cout << "Found map for min but not for max, creating new max map: [" << min << "," << max << "]" << std::endl;
#endif
            auto newRng = std::uniform_int_distribution<>(min, max);
            minMap->second.emplace(max, newRng);
            rng = &newRng;
        }
    }
    else
    {
#if VERBOSE
        std::cout << "Found no map for min and max. Creating new map for both: [" << min << "," << max << "]" << std::endl;
#endif
        auto newRng = std::uniform_int_distribution<>(min, max);
        auto newMap = std::map<int, std::uniform_int_distribution<>>();
        newMap.emplace(max, newRng);
        randomDistributions.emplace(min, newMap);
        rng = &newRng;
    }

    return (*rng)(gen);
}

int rD4() { return rollDice(1, 4); }
int rD6() { return rollDice(1, 6); }
int rD8() { return rollDice(1, 8); }
int rD10() { return rollDice(1, 10); }
int rD12() { return rollDice(1, 12); }
int rD20() { return rollDice(1, 20); }
int rD100() { return rollDice(1, 100); }

struct PointDistribution
{
    std::array<uint8_t, 6> points;
    double mean;
    double variance;

    std::string PointsToString()
    {
        std::stringstream ss;
        ss << "<";
        for (int i = 0; i < 6; ++i)
        {
            if (i != 5)
            {
                ss << std::to_string(points[i]) << ", ";
            }
            else
            {
                ss << std::to_string(points[i]);
            }
        }
        ss << ">";
        return ss.str();
    }

    void ProcessStatistics()
    {
        double sum = std::accumulate(points.begin(), points.end(), 0);
        mean = sum / 6.0;

        // var = sum((X - mean)^2) / (n - 1)
        variance = 0.0;
        for (int k = 0; k < points.size(); ++k)
        {
            variance += std::pow((points[k] - mean), 2);
        }
        variance /= (points.size() - 1.0);
    }
};

struct Measurement
{
    PointDistribution minDist;
    PointDistribution medianDist;
    PointDistribution maxDist;
    std::array<PointDistribution, 3> minDistribution;
    std::array<PointDistribution, 3> per25Distribution;
    std::array<PointDistribution, 3> per50Distribution;
    std::array<PointDistribution, 3> per75Distribution;
    std::array<PointDistribution, 3> maxDistribution;
};

struct HeatMapEntry
{
    uint16_t count = 0;
    PointDistribution dist;
};

struct ProcessedBucket
{
    ProcessedBucket(int sampleCount)
        : sampleCount(sampleCount)
    {
        for (int i = 0; i < overallCountPerPoint.size(); ++i)
        {
            overallCountPerPoint[i] = 0;
        }

        allSamples.resize(sampleCount);
    }

    std::vector<PointDistribution> allSamples;

    std::map<uint32_t, uint16_t> pointDistributionHeatMap; // To get e.g. the 5 most commonly rolled distributions
    std::array<HeatMapEntry, 10> mostCommonlyRolledEntries;
    
    std::array<uint32_t, 21> overallCountPerPoint; // To get a distribution of all point values overall

    Measurement mean;
    Measurement variance;

    int sampleCount;
};

void ProcessIntoMeasurement(ProcessedBucket& bucket, Measurement& measurement, std::function<double(const PointDistribution&)> sortGetter)
{
    std::sort(bucket.allSamples.begin(), bucket.allSamples.end(), [&sortGetter](const PointDistribution& a, const PointDistribution& b) { return sortGetter(a) < sortGetter(b); });
    measurement.minDistribution[0] = bucket.allSamples[0];
    measurement.minDistribution[1] = bucket.allSamples[1];
    measurement.minDistribution[2] = bucket.allSamples[2];

    int perc25 = (int)((double)bucket.sampleCount * 0.25);
    measurement.per25Distribution[0] = bucket.allSamples[perc25 - 1];
    measurement.per25Distribution[1] = bucket.allSamples[perc25 + 0];
    measurement.per25Distribution[2] = bucket.allSamples[perc25 + 1];

    int perc50 = (int)((double)bucket.sampleCount * 0.50);
    measurement.per50Distribution[0] = bucket.allSamples[perc50 - 1];
    measurement.per50Distribution[1] = bucket.allSamples[perc50 + 0];
    measurement.per50Distribution[2] = bucket.allSamples[perc50 + 1];

    int perc75 = (int)((double)bucket.sampleCount * 0.75);
    measurement.per75Distribution[0] = bucket.allSamples[perc75 - 1];
    measurement.per75Distribution[1] = bucket.allSamples[perc75 + 0];
    measurement.per75Distribution[2] = bucket.allSamples[perc75 + 1];

    measurement.maxDistribution[0] = bucket.allSamples[bucket.sampleCount - 3];
    measurement.maxDistribution[1] = bucket.allSamples[bucket.sampleCount - 2];
    measurement.maxDistribution[2] = bucket.allSamples[bucket.sampleCount - 1];

    measurement.minDist = bucket.allSamples[0];
    measurement.medianDist = bucket.allSamples[perc50];
    measurement.maxDist = bucket.allSamples[bucket.sampleCount - 1];
}

void measure(std::stringstream& summaryOutString, int sampleCount, std::function<PointDistribution()> distFunc, std::string name, std::string description)
{
    ProcessedBucket bucket(sampleCount);

    for (int i = 0; i < sampleCount; ++i)
    {
        bucket.allSamples[i] = distFunc();


        std::sort(bucket.allSamples[i].points.begin(), bucket.allSamples[i].points.end(), [](const uint8_t& a, const uint8_t& b) { return a > b; });

        uint32_t key = bucket.allSamples[i].points[0];
        for (int k = 1; k < bucket.allSamples[i].points.size(); ++k)
        {
            key = key << 5;
            key += bucket.allSamples[i].points[k];
        }
        if (bucket.pointDistributionHeatMap.find(key) == bucket.pointDistributionHeatMap.end())
        {
            bucket.pointDistributionHeatMap.emplace(key, 0);
        }
        ++bucket.pointDistributionHeatMap[key];

        for (int k = 0; k < bucket.mostCommonlyRolledEntries.size(); ++k)
        {
            if (bucket.mostCommonlyRolledEntries[k].count < bucket.pointDistributionHeatMap[key])
            {
                bucket.mostCommonlyRolledEntries[k].count = bucket.pointDistributionHeatMap[key];
                bucket.mostCommonlyRolledEntries[k].dist = bucket.allSamples[i];

                if (bucket.allSamples[i].points[0] == bucket.allSamples[i].points[1] == bucket.allSamples[i].points[2] == bucket.allSamples[i].points[3] == bucket.allSamples[i].points[4] == bucket.allSamples[i].points[5])
                {
                    std::cout << "BREAK";
                }

                break;
            }
        }

        for (int k = 0; k < bucket.allSamples[i].points.size(); ++k)
        {
            ++bucket.overallCountPerPoint[bucket.allSamples[i].points[k]];
        }
    }

    ProcessIntoMeasurement(bucket, bucket.mean, [](const PointDistribution& dist) { return dist.mean; });
    ProcessIntoMeasurement(bucket, bucket.variance, [](const PointDistribution& dist) { return dist.variance; });

    std::stringstream outString;
    outString << "# " << name << "\n";
    outString << "\n";
    outString << description;
    outString << "\n";
    outString << "\n";
    outString << "## Means" << "\n";
    outString << "| Type   | Mean | Distribution |" << "\n";
    outString << "| ------ | ---- | ------------ |" << "\n";
    outString << "| Min    | " << bucket.mean.minDist.mean << " | " << bucket.mean.minDist.PointsToString() << " |" << "\n";
    outString << "| Median | " << bucket.mean.medianDist.mean << " | " << bucket.mean.medianDist.PointsToString() << " |" << "\n";
    outString << "| Max    | " << bucket.mean.maxDist.mean << " | " << bucket.mean.maxDist.PointsToString() << " |" << "\n";
    outString << "\n";
    outString << "## Variance" << "\n";
    outString << "| Type   | Variance | Distribution |" << "\n";
    outString << "| ------ | ---- | ------------ |" << "\n";
    outString << "| Min    | " << bucket.variance.minDist.variance << " | " << bucket.variance.minDist.PointsToString() << " |" << "\n";
    outString << "| Median | " << bucket.variance.medianDist.variance << " | " << bucket.variance.medianDist.PointsToString() << " |" << "\n";
    outString << "| Max    | " << bucket.variance.maxDist.variance << " | " << bucket.variance.maxDist.PointsToString() << " |" << "\n";
    outString << "\n";
    outString << "## Most common rolled distributions" << "\n";
    outString << "| Distribution | Rolled n times |" << "\n";
    outString << "| ------------ | -------------- |" << "\n";
    for (int i = 0; i < bucket.mostCommonlyRolledEntries.size(); ++i)
    {
        outString << "| " << bucket.mostCommonlyRolledEntries[i].dist.PointsToString() << " | " << std::to_string(bucket.mostCommonlyRolledEntries[i].count) << " |\n";
    }
    outString << "\n";
    outString << "## Times value was selected" << "\n";
    outString << "Total sample count: " << std::to_string(sampleCount) << "\n";
    outString << "| Value | Selected n times |" << "\n";
    outString << "| ----- | ---------------- |" << "\n";
    for (int i = 1; i < bucket.overallCountPerPoint.size(); ++i)
    {
        outString << "| " << std::to_string(i) << " | " << std::to_string(bucket.overallCountPerPoint[i]) << " |\n";
    }

    auto printMeasurement = [](std::stringstream& outString, Measurement& mes, std::string category) {
            outString << "\n";
            outString << "## Sorted by " << category << "\n";
            outString << "| Category | Example # | Mean | Distribution |" << "\n";
            outString << "| -------- | --------- | ---- | ------------ |" << "\n";
            for (int i = 0; i < mes.minDistribution.size(); ++i)
            {
                outString << "| Minimum | " << std::to_string(i + 1) << " | " << mes.minDistribution[i].mean << " | " << mes.minDistribution[i].PointsToString() << " |\n";
            }
            outString << "| | | | |" << "\n";
            for (int i = 0; i < mes.per25Distribution.size(); ++i)
            {
                outString << "| 25th Percentile | " << std::to_string(i + 1) << " | " << mes.per25Distribution[i].mean << " | " << mes.per25Distribution[i].PointsToString() << " |\n";
            }
            outString << "| | | | |" << "\n";
            for (int i = 0; i < mes.per50Distribution.size(); ++i)
            {
                outString << "| Median | " << std::to_string(i + 1) << " | " << mes.per50Distribution[i].mean << " | " << mes.per50Distribution[i].PointsToString() << " |\n";
            }
            outString << "| | | | |" << "\n";
            for (int i = 0; i < mes.per75Distribution.size(); ++i)
            {
                outString << "| 75th Percentile | " << std::to_string(i + 1) << " | " << mes.per75Distribution[i].mean << " | " << mes.per75Distribution[i].PointsToString() << " |\n";
            }
            outString << "| | | | |" << "\n";
            for (int i = 0; i < mes.maxDistribution.size(); ++i)
            {
                outString << "| Maximum | " << std::to_string(i + 1) << " | " << mes.maxDistribution[i].mean << " | " << mes.maxDistribution[i].PointsToString() << " |\n";
            }
        };

    printMeasurement(outString, bucket.mean, "mean");
    printMeasurement(outString, bucket.variance, "variance");

    summaryOutString << outString.str();

    summaryOutString << "\n$$\\\\[3in]$$\n" << std::endl;

    std::ofstream outfile(name + ".md");
    outfile << outString.str();
    outfile.flush();
    outfile.close();

    std::cout << outString.str();
    std::cout << "\n\n\n\n";
}

PointDistribution rD20DistFunc()
{
    PointDistribution sample;
    for (int k = 0; k < sample.points.size(); ++k)
    {
        uint8_t point = rD20();
        sample.points[k] = point;
    }
    sample.ProcessStatistics();

    return sample;
}

PointDistribution standardPoints()
{
    PointDistribution sample;
    sample.points[0] = 15;
    sample.points[1] = 14;
    sample.points[2] = 12;
    sample.points[3] = 11;
    sample.points[4] = 10;
    sample.points[5] = 8;
    sample.ProcessStatistics();
    return sample;
}

PointDistribution gridAllDirs()
{
    PointDistribution sample;
    std::array<std::array<uint8_t, 3>, 3> grid;
    for (int i = 0; i < grid.size(); ++i)
    {
        for (int k = 0; k < grid[0].size(); ++k)
        {
            grid[i][k] = rD6();
        }
    }

    std::array<uint8_t, 8> numbers;
    // horizontal
    numbers[0] = grid[0][0] + grid[1][0] + grid[2][0];
    numbers[1] = grid[0][1] + grid[1][1] + grid[2][1];
    numbers[2] = grid[0][2] + grid[1][2] + grid[2][2];

    // vertical
    numbers[3] = grid[0][0] + grid[0][1] + grid[0][2];
    numbers[4] = grid[1][0] + grid[1][1] + grid[1][2];
    numbers[5] = grid[2][0] + grid[2][1] + grid[2][2];

    // diagonal
    numbers[6] = grid[0][0] + grid[1][1] + grid[2][2];
    numbers[7] = grid[0][2] + grid[1][1] + grid[2][0];

    std::sort(numbers.begin(), numbers.end());

    for (int i = 0; i < sample.points.size(); ++i)
    {
        sample.points[i] = numbers[numbers.size() - 1 - i];
    }
    sample.ProcessStatistics();
    return sample;
}

PointDistribution gridAllDirsCenter3()
{
    PointDistribution sample;
    std::array<std::array<uint8_t, 3>, 3> grid;
    for (int i = 0; i < grid.size(); ++i)
    {
        for (int k = 0; k < grid[0].size(); ++k)
        {
            grid[i][k] = rD6();
        }
    }

    grid[1][1] = 3;

    std::array<uint8_t, 8> numbers;
    // horizontal
    numbers[0] = grid[0][0] + grid[1][0] + grid[2][0];
    numbers[1] = grid[0][1] + grid[1][1] + grid[2][1];
    numbers[2] = grid[0][2] + grid[1][2] + grid[2][2];

    // vertical
    numbers[3] = grid[0][0] + grid[0][1] + grid[0][2];
    numbers[4] = grid[1][0] + grid[1][1] + grid[1][2];
    numbers[5] = grid[2][0] + grid[2][1] + grid[2][2];

    // diagonal
    numbers[6] = grid[0][0] + grid[1][1] + grid[2][2];
    numbers[7] = grid[0][2] + grid[1][1] + grid[2][0];

    std::sort(numbers.begin(), numbers.end());

    for (int i = 0; i < sample.points.size(); ++i)
    {
        sample.points[i] = numbers[numbers.size() - 1 - i];
    }
    sample.ProcessStatistics();
    return sample;
}

PointDistribution gridAllDirsFlip()
{
    PointDistribution sample;
    std::array<std::array<uint8_t, 3>, 3> grid;
    for (int i = 0; i < grid.size(); ++i)
    {
        for (int k = 0; k < grid[0].size(); ++k)
        {
            grid[i][k] = rD6();
        }
    }

    for (int i = 0; i < sample.points.size(); ++i)
    {
        std::array<uint8_t, 8> numbers;
        // horizontal
        numbers[0] = grid[0][0] + grid[1][0] + grid[2][0];
        numbers[1] = grid[0][1] + grid[1][1] + grid[2][1];
        numbers[2] = grid[0][2] + grid[1][2] + grid[2][2];

        // vertical
        numbers[3] = grid[0][0] + grid[0][1] + grid[0][2];
        numbers[4] = grid[1][0] + grid[1][1] + grid[1][2];
        numbers[5] = grid[2][0] + grid[2][1] + grid[2][2];

        // diagonal
        numbers[6] = grid[0][0] + grid[1][1] + grid[2][2];
        numbers[7] = grid[0][2] + grid[1][1] + grid[2][0];

        uint8_t maxVal = 0;
        std::array<uint8_t, numbers.size()> sectionsWithMaxValue;
        uint8_t equalSectionCount = 0;

        for (int i = 0; i < numbers.size(); ++i)
        {
            if (numbers[i] > maxVal)
            {
                equalSectionCount = 1;
                sectionsWithMaxValue[0] = i;
                maxVal = numbers[i];
            }
            else if (numbers[i] == maxVal)
            {
                ++equalSectionCount;
                sectionsWithMaxValue[equalSectionCount - 1] = i;
            }
        }

        sample.points[i] = maxVal;

        uint8_t chosenSection = 0;
        if (equalSectionCount > 1)
        {
            chosenSection = rollDice(0, equalSectionCount - 1);
        }

        auto flip = [](uint8_t& num1, uint8_t& num2, uint8_t& num3) {
            num1 = 7 - num1;
            num2 = 7 - num2;
            num3 = 7 - num3;
        };

        switch (sectionsWithMaxValue[chosenSection])
        {
            case 0:
                flip(grid[0][0], grid[1][0], grid[2][0]);
                break;
            case 1:
                flip(grid[0][1], grid[1][1], grid[2][1]);
                break;
            case 2:
                flip(grid[0][2], grid[1][2], grid[2][2]);
                break;
            case 3:
                flip(grid[0][0], grid[0][1], grid[0][2]);
                break;
            case 4:
                flip(grid[1][0], grid[1][1], grid[1][2]);
                break;
            case 5:
                flip(grid[2][0], grid[2][1], grid[2][2]);
                break;
            case 6:
                flip(grid[0][0], grid[1][1], grid[2][2]);
                break;
            case 7:
                flip(grid[0][2], grid[1][1], grid[2][0]);
                break;
        }
    }

    sample.ProcessStatistics();
    return sample;
}

PointDistribution gridNoDiagonalsFlip()
{
    PointDistribution sample;
    std::array<std::array<uint8_t, 3>, 3> grid;
    for (int i = 0; i < grid.size(); ++i)
    {
        for (int k = 0; k < grid[0].size(); ++k)
        {
            grid[i][k] = rD6();
        }
    }

    for (int i = 0; i < sample.points.size(); ++i)
    {
        std::array<uint8_t, 6> numbers;
        // horizontal
        numbers[0] = grid[0][0] + grid[1][0] + grid[2][0];
        numbers[1] = grid[0][1] + grid[1][1] + grid[2][1];
        numbers[2] = grid[0][2] + grid[1][2] + grid[2][2];

        // vertical
        numbers[3] = grid[0][0] + grid[0][1] + grid[0][2];
        numbers[4] = grid[1][0] + grid[1][1] + grid[1][2];
        numbers[5] = grid[2][0] + grid[2][1] + grid[2][2];

        uint8_t maxVal = 0;
        std::array<uint8_t, numbers.size()> sectionsWithMaxValue;
        uint8_t equalSectionCount = 0;

        for (int i = 0; i < numbers.size(); ++i)
        {
            if (numbers[i] > maxVal)
            {
                equalSectionCount = 1;
                sectionsWithMaxValue[0] = i;
                maxVal = numbers[i];
            }
            else if (numbers[i] == maxVal)
            {
                ++equalSectionCount;
                sectionsWithMaxValue[equalSectionCount - 1] = i;
            }
        }

        sample.points[i] = maxVal;

        uint8_t chosenSection = 0;
        if (equalSectionCount > 1)
        {
            chosenSection = rollDice(0, equalSectionCount - 1);
        }

        auto flip = [](uint8_t& num1, uint8_t& num2, uint8_t& num3) {
            num1 = 7 - num1;
            num2 = 7 - num2;
            num3 = 7 - num3;
            };

        switch (sectionsWithMaxValue[chosenSection])
        {
        case 0:
            flip(grid[0][0], grid[1][0], grid[2][0]);
            break;
        case 1:
            flip(grid[0][1], grid[1][1], grid[2][1]);
            break;
        case 2:
            flip(grid[0][2], grid[1][2], grid[2][2]);
            break;
        case 3:
            flip(grid[0][0], grid[0][1], grid[0][2]);
            break;
        case 4:
            flip(grid[1][0], grid[1][1], grid[1][2]);
            break;
        case 5:
            flip(grid[2][0], grid[2][1], grid[2][2]);
            break;
        }
    }

    sample.ProcessStatistics();
    return sample;
}

PointDistribution r4D6D1()
{
    PointDistribution sample;

    auto rollOnce = []() {
        std::array<uint8_t, 4> rolls;
        for (int i = 0; i < rolls.size(); ++i)
        {
            rolls[i] = rD6();
        }

        std::sort(rolls.begin(), rolls.end());
        return std::accumulate(rolls.begin() + 1, rolls.end(), 0);
    };

    for (int i = 0; i < sample.points.size(); ++i)
    {
        sample.points[i] = rollOnce();
    }

    sample.ProcessStatistics();
    return sample;
}

int main()
{
    std::stringstream summaryOutString;
    //measure(summaryOutString, 10000, rD20DistFunc, "D20 Raw", "Just roll a d20 and take 6 numbers.");
    //measure(summaryOutString, 10, standardPoints, "Standard Points", "Standard distribution of <15 14 12 11 10 8>.");
    //measure(summaryOutString, 10000, gridAllDirs, "Grid All Directions", "Roll 9 d6. Lay in grid 3 x 3. Sum all horizontal, vertical and diagonal lines. Take the highest 6.");
    //measure(summaryOutString, 10000, gridAllDirsCenter3, "Grid All Directions Center 3", "Roll 8 d6. Lay in grid 3 x 3 and add a 3 as the center value. Sum all horizontal, vertical and diagonal lines. Take the highest 6.");
    //measure(summaryOutString, 10000, gridAllDirsFlip, "Grid All Directions Flip", "Roll 9 d6. Lay in grid 3 x 3. Sum all horizontal, vertical and diagonal lines. Take the highest. Flip the dice of the chosen lane. Repeat 6 times.");
    //measure(summaryOutString, 10000, gridNoDiagonalsFlip, "Grid No Diagonals Flip", "Roll 9 d6. Lay in grid 3 x 3. Sum all horizontal and vertical lines. Take the highest. Flip the dice of the chosen lane. Repeat 6 times.");
    measure(summaryOutString, 10000, r4D6D1, "Roll 4D6D1", "Roll 4 d6. Drop the lowest. Take sum. Repeat 6 times.");

    std::ofstream outfile("allMethods.md");
    outfile << summaryOutString.str();
    outfile.flush();
    outfile.close();
}