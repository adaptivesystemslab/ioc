function model = inputData(model, q, dq, ddq)
    model.position = q;
    model.velocity = dq;
    model.acceleration = ddq;
end