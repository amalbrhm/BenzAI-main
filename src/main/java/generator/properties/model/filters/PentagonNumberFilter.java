package generator.properties.model.filters;

import benzenoid.Benzenoid;
import benzenoid.CycleType;
import generator.properties.model.ModelPropertySet;
import generator.properties.model.expression.BinaryNumericalExpression;
import generator.properties.model.expression.PropertyExpression;

import java.util.ArrayList;

public class PentagonNumberFilter extends Filter {

    @Override
    public boolean test(Benzenoid molecule, ArrayList<PropertyExpression> propertyExpressionList, ModelPropertySet modelPropertySet) {
        int pentagonCount = 0;
        // TODO: Implement robust pentagon counting in Benzenoid class or here.
        // Using CycleType as a placeholder, assuming getNbHexagons() might give total face count
        // and getHexagonCycleType(i) can identify a C5 face.
        // This needs verification based on how Benzenoid handles non-hexagonal structures.
        if (molecule != null) { // Ensure molecule is not null
            for (int i = 0; i < molecule.getNbHexagons(); i++) {
                if (molecule.getHexagonCycleType(i) == CycleType.C5) {
                    pentagonCount++;
                }
            }
        }


        for (PropertyExpression expression : propertyExpressionList) {
            if (expression instanceof BinaryNumericalExpression) {
                BinaryNumericalExpression numExpr = (BinaryNumericalExpression) expression;
                if (!numExpr.test(pentagonCount, numExpr.getOperator(), numExpr.getValue())) {
                    return false;
                }
            }
        }
        return true;
    }
}