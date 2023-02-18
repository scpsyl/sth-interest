#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define INPUT_SIZE 1  // 输入层大小
#define HIDDEN_SIZE 3 // 隐藏层大小
#define OUTPUT_SIZE 1 // 输出层大小

#define LEARNING_RATE 0.1 // 学习率
//#define NUM_ITERATIONS_1  10 // 迭代次数
#define NUM_ITERATIONS_2  100
#define MOMENTUM 0.9      // 动量参数

#define INTEGRAL_MAX 1000 // 积分限幅
#define INTEGRAL_MIN -1000

typedef struct
{
    double *inputs;  // 输入
    double *weights; // 权重
    double bias;     // 偏差
    double output;   // 输出
} Neuron;

typedef struct
{
    Neuron *neurons; // 神经元
    int size;        // 大小
} Layer;

typedef struct
{
    Layer inputLayer;  // 输入层
    Layer hiddenLayer; // 隐藏层
    Layer outputLayer; // 输出层
} NeuralNetwork;

typedef struct
{
    double Kp;            // 比例常数
    double Ki;            // 积分常数
    double Kd;            // 微分常数
    double integral;      // 积分
    double previousError; // 上一次误差
} PIDController;

void initializeNeuron(Neuron *neuron, int inputSize)
{
    neuron->inputs = (double *)malloc(sizeof(double) * inputSize);
    neuron->weights = (double *)malloc(sizeof(double) * inputSize);
    neuron->bias = (double)rand() / RAND_MAX * 2.0 - 1.0;
    neuron->output = 0.0;
    for (int i = 0; i < inputSize; i++)
    {
        neuron->weights[i] = (double)rand() / RAND_MAX * 2.0 - 1.0;
    }
}

void initializeLayer(Layer *layer, int size, int inputSize)
{
    layer->neurons = (Neuron *)malloc(sizeof(Neuron) * size);
    layer->size = size;
    for (int i = 0; i < size; i++)
    {
        initializeNeuron(&layer->neurons[i], inputSize);
    }
}

void initializeNetwork(NeuralNetwork *network)
{
    initializeLayer(&network->inputLayer, INPUT_SIZE, 0);
    initializeLayer(&network->hiddenLayer, HIDDEN_SIZE, INPUT_SIZE);
    initializeLayer(&network->outputLayer, OUTPUT_SIZE, HIDDEN_SIZE);
}

double sigmoid(double x)
{
    return 1.0 / (1.0 + exp(-x));
}

void forward(NeuralNetwork *network, double *inputs)
{
    // 将输入层输出设置为输入
    for (int i = 0; i < INPUT_SIZE; i++)
    {
        network->inputLayer.neurons[i].output = inputs[i];
    }

    // 计算隐藏层输出
    for (int i = 0; i < HIDDEN_SIZE; i++)
    {
        double sum = 0.0;
        for (int j = 0; j < INPUT_SIZE; j++)
        {
            sum += network->inputLayer.neurons[j].output * network->hiddenLayer.neurons[i].weights[j];
        }
        sum += network->hiddenLayer.neurons[i].bias;
        network->hiddenLayer.neurons[i].output = sigmoid(sum);
    }
}
double calculateError(double target, double output)
{
    return target - output;
}

double calculateDerivative(double output)
{
    return output * (1.0 - output);
}
void backward(NeuralNetwork *network, double *targets)
{
    // 计算输出层误差
    for (int i = 0; i < OUTPUT_SIZE; i++)
    {
        Neuron *outputNeuron = &network->outputLayer.neurons[i];
        double error = calculateError(targets[i], outputNeuron->output);
        double derivative = calculateDerivative(outputNeuron->output);
        outputNeuron->output = error * derivative;
    }
    // 计算隐藏层误差
    for (int i = 0; i < HIDDEN_SIZE; i++)
    {
        Neuron *hiddenNeuron = &network->hiddenLayer.neurons[i];
        double sum = 0.0;
        for (int j = 0; j < OUTPUT_SIZE; j++)
        {
            sum += network->outputLayer.neurons[j].output * network->outputLayer.neurons[j].weights[i];
        }
        double derivative = calculateDerivative(hiddenNeuron->output);
        hiddenNeuron->output = sum * derivative;
    }

    // 更新输出层权重和偏差
    for (int i = 0; i < OUTPUT_SIZE; i++)
    {
        Neuron *outputNeuron = &network->outputLayer.neurons[i];
        for (int j = 0; j < HIDDEN_SIZE; j++)
        {
            double delta = LEARNING_RATE * outputNeuron->output * network->hiddenLayer.neurons[j].output;
            outputNeuron->weights[j] += delta;
        }
        outputNeuron->bias += LEARNING_RATE * outputNeuron->output;
    }

    // 更新隐藏层权重和偏差
    for (int i = 0; i < HIDDEN_SIZE; i++)
    {
        Neuron *hiddenNeuron = &network->hiddenLayer.neurons[i];
        for (int j = 0; j < INPUT_SIZE; j++)
        {
            double delta = LEARNING_RATE * hiddenNeuron->output * network->inputLayer.neurons[j].output;
            hiddenNeuron->weights[j] += delta;
        }
        hiddenNeuron->bias += LEARNING_RATE * hiddenNeuron->output;
    }
}

void train(NeuralNetwork *network, double *inputs, double *targets)
{
    forward(network, inputs);
    backward(network, targets);
}
void initializePIDController(PIDController *pid, double Kp, double Ki, double Kd)
{
    pid->Kp = Kp;
    pid->Ki = Ki;
    pid->Kd = Kd;
    pid->integral = 0.0;
    pid->previousError = 0.0;
}

double calculatePIDOutput(PIDController *pid, double error)
{
    double proportional = pid->Kp * error;
    pid->integral += error;
    double integral = pid->Ki * pid->integral;

    // 积分限制
    if (pid->integral > INTEGRAL_MAX)
    {
        pid->integral = INTEGRAL_MAX;
    }
    else if (pid->integral < INTEGRAL_MIN)
    {
        pid->integral = INTEGRAL_MIN;
    }
    double derivative = pid->Kd * (error - pid->previousError);
    pid->previousError = error;

    double output = proportional + pid->integral + derivative;
    return output;
}

int main()
{
    NeuralNetwork network;
    initializeNetwork(&network);
    PIDController pid;
    initializePIDController(&pid, 0.1, 0.01, 0.01);

    double setpoint = 30.0;
    double processVariable = 0.0;
#ifdef NUM_ITERATIONS_1
    for (int i = 0; i < NUM_ITERATIONS_1; i++)
#endif

#ifdef NUM_ITERATIONS_2
    for (int i = 0; i < NUM_ITERATIONS_2; i++)
#endif
    {
        double error = setpoint - processVariable;
        double output = calculatePIDOutput(&pid, error);
        processVariable += output;

        double inputs[INPUT_SIZE] = {processVariable};
        double targets[OUTPUT_SIZE] = {setpoint};

        train(&network, inputs, targets);
    }

    printf("Final output: %f\n", processVariable);
    printf("Loss: %f\n",(setpoint-processVariable));

    return 0;
}
