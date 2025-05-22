package generator.properties.model.filters;

import benzenoid.Benzenoid;
import benzenoid.CycleType;
import generator.properties.model.ModelPropertySet;
import generator.properties.model.expression.BinaryNumericalExpression;
import generator.properties.model.expression.PropertyExpression;

import java.util.ArrayList;

public class HeptagonNumberFilter extends Filter {

    @Override
    public boolean test(Benzenoid molecule, ArrayList<PropertyExpression> propertyExpressionList, ModelPropertySet modelPropertySet) {
        int heptagonCount = 0;
        // TODO: Implement robust heptagon counting in Benzenoid class or here.
        if (molecule != null) { // Ensure molecule is not null
            for (int i = 0; i < molecule.getNbHexagons(); i++) {
                if (molecule.getHexagonCycleType(i) == CycleType.C7) {
                    heptagonCount++;
                }
            }
        }

        for (PropertyExpression expression : propertyExpressionList) {
            if (expression instanceof BinaryNumericalExpression) {
                BinaryNumericalExpression numExpr = (BinaryNumericalExpression) expression;
                if (!numExpr.test(heptagonCount, numExpr.getOperator(), numExpr.getValue())) {
                    return false;
                }
            }
        }
        return true;
    }
}