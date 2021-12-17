function [accuracy, sensitivity, specificity, precision, dice, loU] = criterion(predict, label)
TP = sum(predict(:) ~= 0 & label(:) ~= 0);
TN = sum(predict(:) == 0 & label(:) == 0);
FP = sum(predict(:) ~= 0 & label(:) == 0);
FN = sum(predict(:) == 0 & label(:) ~= 0);
accuracy = (TP + TN) / numel(label); % 精度
sensitivity = TP / (TP + FN); % 敏感性/召回率recall
specificity = TN / (TN + FP); % 特异性
precision = TP / (TP + FP); % 查准率：(预测为1且正确预测的样本数)/(所有预测为1的样本数)
dice = 2 * TP / (2 * TP + FN + FP);
loU = TP / (FP + TP + FN); % ???????????? ???? ?????
fprintf('accuracy:%g\tsensitivity:%g\tspecificity:%g\tprecision:%g\tdice:%g\tloU:%g\n', ...
    accuracy, sensitivity, specificity, precision, dice, loU);
end