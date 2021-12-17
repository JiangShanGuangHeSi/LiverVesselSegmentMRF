function [accuracy, sensitivity, specificity, precision, dice, loU] = criterion(predict, label)
TP = sum(predict(:) ~= 0 & label(:) ~= 0);
TN = sum(predict(:) == 0 & label(:) == 0);
FP = sum(predict(:) ~= 0 & label(:) == 0);
FN = sum(predict(:) == 0 & label(:) ~= 0);
accuracy = (TP + TN) / numel(label); % ����
sensitivity = TP / (TP + FN); % ������/�ٻ���recall
specificity = TN / (TN + FP); % ������
precision = TP / (TP + FP); % ��׼�ʣ�(Ԥ��Ϊ1����ȷԤ���������)/(����Ԥ��Ϊ1��������)
dice = 2 * TP / (2 * TP + FN + FP);
loU = TP / (FP + TP + FN); % ???????????? ???? ?????
fprintf('accuracy:%g\tsensitivity:%g\tspecificity:%g\tprecision:%g\tdice:%g\tloU:%g\n', ...
    accuracy, sensitivity, specificity, precision, dice, loU);
end