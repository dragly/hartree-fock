class MockBasisFunction : public BasisFunction {
 public:
  MOCK_METHOD1(importantFunction,
      void(int someNumber));
};
